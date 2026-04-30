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
        help="Orchestrate coverage.py analysis during e2e tests through a special Docker Compose overlay that runs it in the FastAPI and Celery worker containers.",
    )
