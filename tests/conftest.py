"""
Top-level pytest configuration for DivBase
"""

from structlog.typing import EventDict

REGRESSION_GUARD_PREFIX = "Regression guard failed:"


def pytest_addoption(parser):
    """Custom command-line options to pytest. In this case to run tests marked as slow"""
    parser.addoption(
        "--run-slow",
        action="store_true",
        default=False,
        help="Run slow tests",
    )
    parser.addoption(
        "--coverage-docker",
        action="store_true",
        default=False,
        help="Enable running test coverage analysis for all tests, including those that run code inside docker containers. Do not use this flag directly, use the `scripts/run_tests_with_coverage.sh` to do coverage analysis",
    )


def _text_in_logs(text: str, logs: list[EventDict]) -> bool:
    """
    Helper fn to check if a specific chunk of text is present in any of the log msgs/events.

    This works with structlog's EventDict log structure.
    Used in both e2e and unit tests, hence why it is here.
    """
    return any(text in log.get("event", "") for log in logs)
