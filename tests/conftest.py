"""
Top-level pytest configuration for DivBase
"""

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
        "--run-very-slow",
        action="store_true",
        default=False,
        help="Run tests marked very-slow (may take several minutes). Also runs tests marked slow.",
    )
    parser.addoption(
        "--coverage-docker",
        action="store_true",
        default=False,
        help="Enable running test coverage analysis for all tests, including those that run code inside docker containers. Do not use this flag directly, use the `scripts/run_tests_with_coverage.sh` to do coverage analysis",
    )


def pytest_configure(config):
    """When --run-very-slow is set, also enable --run-slow."""
    if config.getoption("--run-very-slow", default=False):
        config.option.run_slow = True
