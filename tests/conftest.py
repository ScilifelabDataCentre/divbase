"""
Top-level pytest configuration for the DivBase project.
It handles spinning up the job system docker stack for the duration of the test session, and the tear-down afterwards.
"""

import os

import pytest

from tests.helpers.docker_testing_stack_setup import start_compose_stack, stop_compose_stack
from tests.helpers.minio_setup import MINIO_FAKE_ACCESS_KEY, MINIO_FAKE_SECRET_KEY, setup_minio_data


@pytest.fixture(autouse=True, scope="session")
def set_env_vars():
    """
    Set environment variables for duration of the entire test session.

    For Celery broker and result backend to run in the
    testing docker stack defined in job_system_compose_tests.yaml.
    Separated from job_system_docker_stack fixture to ensure that these env variables are
    used in all tests that require Celery, even if the job system stack was started outside of pytest.
    """
    os.environ["CELERY_BROKER_URL"] = "pyamqp://guest@localhost:5673//"
    os.environ["CELERY_RESULT_BACKEND"] = "redis://localhost:6380/0"
    os.environ["FLOWER_USER"] = "floweradmin"
    os.environ["FLOWER_PASSWORD"] = "badpassword"
    os.environ["FLOWER_BASE_URL"] = "http://localhost:5556"

    os.environ["DIVBASE_S3_ACCESS_KEY"] = MINIO_FAKE_ACCESS_KEY
    os.environ["DIVBASE_S3_SECRET_KEY"] = MINIO_FAKE_SECRET_KEY


@pytest.fixture(autouse=True, scope="session")
def docker_testing_stack():
    """
    Start job system docker stack, and stop after all tests run.
    """
    try:
        start_compose_stack()
        setup_minio_data()
        yield
    finally:
        stop_compose_stack()
