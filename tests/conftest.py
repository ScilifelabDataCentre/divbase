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
        "--run-very-slow",
        action="store_true",
        default=False,
        help="Run tests marked very-slow (may take several minutes). Also runs tests marked slow.",
    )


def pytest_configure(config):
    """When --run-very-slow is set, also enable --run-slow."""
    if config.getoption("--run-very-slow", default=False):
        config.option.run_slow = True
