"""
Unit tests for the uvicorn access log filter that drops healthy healthcheck requests.
"""

import logging

import pytest

from divbase_api.logging_config import SkipHealthyApiChecksFilter


def _access_record(path: str, status: int, client: str = "1.2.3.4:5000") -> logging.LogRecord:
    """Build a log record shaped like the ones uvicorn emits on the uvicorn.access logger."""
    return logging.LogRecord(
        name="uvicorn.access",
        level=logging.INFO,
        pathname=__file__,
        lineno=0,
        msg='%s - "%s %s HTTP/%s" %d',
        args=(client, "GET", path, "1.1", status),
        exc_info=None,
    )


@pytest.mark.parametrize(
    "client, path, status, kept",
    [
        ("1.2.3.4:5000", "/api/v1/core/health", 200, False),
        ("1.2.3.4:5000", "/api/v1/core/health?x=1", 200, False),
        # "200" elsewhere in the line must not hide a failed check
        ("10.0.0.5:52003", "/api/v1/core/health", 503, True),
        ("10.200.0.1:1", "/api/v1/core/health", 500, True),
        ("1.2.3.4:5000", "/api/v1/core/health?x=200", 404, True),
        ("1.2.3.4:5000", "/api/v1/core/healthz", 200, True),
        ("1.2.3.4:5000", "/api/v1/other", 200, True),
    ],
)
def test_health_filter_decision(client, path, status, kept):
    record = _access_record(path, status, client)
    assert SkipHealthyApiChecksFilter().filter(record) is kept
