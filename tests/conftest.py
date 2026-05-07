"""
Top-level pytest configuration for DivBase
"""


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
