"""
Structlog configuration for the divbase-api.

- local-dev/test environments: ConsoleRenderer
- All other environments: JSON for Loki.

Usage in all modules is simply:
logger = structlog.get_logger(__name__)
(but if standard logging is used in some places, it will still work fine)
"""

import logging
import sys

import structlog

from divbase_lib.divbase_constants import LOCAL_DEV_ENVIRONMENTS


def configure_logging(log_level: str, environment: str) -> None:
    """
    At app startup, configure structlogging

    This covers logs from both:
    1. structlog logs (logs from calls to structlog.get_logger())
    2. non-structlog logs (e.g. 3rd party lib logs like uvicorn, sqlalchemy, or if there are standard logging calls from our codebase)

    See this guide:
    https://wazaari.dev/blog/fastapi-structlog-integration
    """
    if environment in LOCAL_DEV_ENVIRONMENTS:
        renderer: structlog.types.Processor = structlog.dev.ConsoleRenderer(colors=False)
    else:
        renderer = structlog.processors.JSONRenderer()

    level = getattr(logging, log_level.upper())

    shared_processors: list[structlog.types.Processor] = [
        structlog.contextvars.merge_contextvars,
        structlog.stdlib.add_logger_name,
        structlog.stdlib.add_log_level,
        structlog.stdlib.PositionalArgumentsFormatter(),
        structlog.processors.TimeStamper(fmt="iso"),
        structlog.processors.StackInfoRenderer(),
    ]

    structlog.configure(
        processors=[
            *shared_processors,
            structlog.stdlib.ProcessorFormatter.wrap_for_formatter,
        ],
        wrapper_class=structlog.make_filtering_bound_logger(level),
        context_class=dict,
        logger_factory=structlog.stdlib.LoggerFactory(),
        cache_logger_on_first_use=True,
    )

    # handles non structlog logs
    formatter = structlog.stdlib.ProcessorFormatter(
        foreign_pre_chain=shared_processors,
        processors=[
            structlog.stdlib.ProcessorFormatter.remove_processors_meta,
            renderer,
        ],
    )

    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(formatter)

    root_logger = logging.getLogger()
    root_logger.handlers = [handler]
    root_logger.setLevel(level)

    # prevents double logging by uvicorn
    for name in ("uvicorn", "uvicorn.error", "uvicorn.access", "uvicorn.asgi"):
        uvicorn_logger = logging.getLogger(name)
        uvicorn_logger.handlers.clear()
        uvicorn_logger.propagate = True
