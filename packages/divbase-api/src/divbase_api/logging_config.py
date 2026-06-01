"""
Structlog configuration for the divbase-api.

- local-dev/test environments: ConsoleRenderer
- All other environments: JSON for Loki.

Usage in all modules is simply:
logger = structlog.get_logger(__name__)
(but if standard logging is used in some places, it will still work fine)
"""

import logging
import socket
import sys
from logging.handlers import RotatingFileHandler
from pathlib import Path

import structlog
from structlog.types import EventDict

from divbase_lib.divbase_constants import LOCAL_DEV_ENVIRONMENTS

LOG_FILES_DIR = Path("/logs")

# shared for both (1) application logs and (2) task logs given to users
_SHARED_PROCESSORS: list[structlog.types.Processor] = [
    structlog.contextvars.merge_contextvars,
    structlog.stdlib.add_log_level,
    structlog.stdlib.add_logger_name,
    structlog.stdlib.PositionalArgumentsFormatter(),
    structlog.processors.StackInfoRenderer(),
]

_APP_LOG_PROCESSORS: list[structlog.types.Processor] = [
    *_SHARED_PROCESSORS,
    structlog.processors.TimeStamper(fmt="iso"),
]


class SkipHealthyApiChecksFilter(logging.Filter):
    """
    Filter out uvicorn access logs for the healthcheck endpoint from being stored/logged (otherwise noisy)
    only those that are healthy are filtered, and this is not applied to the stream handler, just the file handler (if used).
    """

    def filter(self, record: logging.LogRecord) -> bool:
        if record.name != "uvicorn.access":
            return True

        message = record.getMessage()
        return "/api/v1/core/health" not in message or "200" not in message


def configure_logging(log_level: str, environment: str, service_name: str, log_to_file: bool = False) -> None:
    """
    At app startup, configure structlogging

    This covers logs from both:
    1. structlog logs (logs from calls to structlog.get_logger())
    2. non-structlog logs (e.g. 3rd party lib logs like uvicorn, sqlalchemy, or if there are standard logging calls from our codebase)

    You can optionally also log to file (in addition to stdout).

    See this guide:
    https://wazaari.dev/blog/fastapi-structlog-integration
    """
    if environment in LOCAL_DEV_ENVIRONMENTS:
        renderer: structlog.types.Processor = structlog.dev.ConsoleRenderer(
            colors=False, exception_formatter=structlog.dev.plain_traceback
        )
    else:
        renderer = structlog.processors.JSONRenderer()

    level = getattr(logging, log_level.upper())

    structlog.configure(
        processors=[
            *_APP_LOG_PROCESSORS,
            structlog.stdlib.ProcessorFormatter.wrap_for_formatter,
        ],
        wrapper_class=structlog.make_filtering_bound_logger(level),
        context_class=dict,
        logger_factory=structlog.stdlib.LoggerFactory(),
        cache_logger_on_first_use=True,
    )

    # handles non structlog logs
    formatter = structlog.stdlib.ProcessorFormatter(
        foreign_pre_chain=_APP_LOG_PROCESSORS,
        processors=[
            structlog.stdlib.ProcessorFormatter.remove_processors_meta,
            renderer,
        ],
    )

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)

    handlers: list[logging.Handler] = [stream_handler]
    if log_to_file:
        instance_name = socket.gethostname()
        log_file = LOG_FILES_DIR / f"{service_name}-{instance_name}.log"

        file_handler = RotatingFileHandler(
            filename=log_file,
            maxBytes=10 * 1024 * 1024,  # 10 MB
            backupCount=2,
            encoding="utf-8",
        )
        file_handler.setFormatter(formatter)
        if service_name == "divbase_api":
            # don't store 200 health check logs in the file
            file_handler.addFilter(SkipHealthyApiChecksFilter())

        handlers.append(file_handler)

    root_logger = logging.getLogger()
    root_logger.handlers = handlers
    root_logger.setLevel(level)

    # prevents double logging by uvicorn
    for name in ("uvicorn", "uvicorn.error", "uvicorn.access", "uvicorn.asgi"):
        uvicorn_logger = logging.getLogger(name)
        uvicorn_logger.handlers.clear()
        uvicorn_logger.propagate = True


def create_user_task_log_handler(log_file: Path) -> logging.Handler:
    """
    Create an empty log file for user query logs and return a FileHandler that can writes human-readable log lines to it.
    Note that as these are for user facing logs we:
        1. drop stacktraces (they could leak something sensitive) and
        2. simplify log lines by removing some fields like date, module used etc...
    """
    log_file.parent.mkdir(exist_ok=True, parents=True)
    formatter = structlog.stdlib.ProcessorFormatter(
        foreign_pre_chain=_SHARED_PROCESSORS,
        processors=[
            structlog.stdlib.ProcessorFormatter.remove_processors_meta,
            _simplify_user_log_events,
            _strip_exc_info,
            structlog.dev.ConsoleRenderer(colors=False),
        ],
    )
    handler = logging.FileHandler(log_file, mode="w", encoding="utf-8")
    handler.setFormatter(formatter)
    return handler


def _simplify_user_log_events(_logger, _method, event_dict: EventDict) -> EventDict:
    """
    Remove what would be a lot of noise from user-facing task log files.
    It's same task_id + name throughout a user's log file and timestarted is shown in the log header.
    'logger' is the modules name, decided to also remove that for simplicity
    """
    event_dict.pop("timestamp", None)
    event_dict.pop("task_id", None)
    event_dict.pop("task_name", None)
    event_dict.pop("logger", None)
    return event_dict


def _strip_exc_info(_logger, _method, event_dict: EventDict) -> EventDict:
    """Strip traceback info from user-facing log lines — show the error message only."""
    event_dict.pop("exc_info", None)
    event_dict.pop("exception", None)
    return event_dict
