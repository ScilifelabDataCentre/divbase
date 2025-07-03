"""
Top-level pytest configuration for the DivBase project.
It handles spinning up the job system docker stack for the duration of the test session, and the tear-down afterwards.
"""

import os

import pytest

from tests.helpers.job_system_setup import start_job_system, stop_job_system
from tests.helpers.minio_setup import (
    setup_minio_data,
    start_minio,
    stop_minio,
)


@pytest.fixture(autouse=True, scope="session")
def set_celery_env_vars():
    """
    Set environment variables for Celery broker and result backend to run in the
    testing docker stack defined in job_system_compose_tests.yaml.
    Separated from job_system_docker_stack fixture to ensure that these env variables are
    used in all tests that require Celery, even if the job system stack was started outside of pytest.
    """
    os.environ["CELERY_BROKER_URL"] = "pyamqp://guest@localhost:5673//"
    os.environ["CELERY_RESULT_BACKEND"] = "redis://localhost:6380/0"


@pytest.fixture(autouse=True, scope="session")
def job_system_docker_stack():
    """
    Start job system docker stack, and stop after all tests run.
    (The print message is only displayed when pytest is run with the -s option)
    """
    try:
        start_job_system()
        yield
    finally:
        stop_job_system()


@pytest.fixture(scope="session")
def minio_server():
    """Start Minio server, set up test buckets, users etc... and stop after all tests run"""
    try:
        start_minio()
        setup_minio_data()
        yield "http://localhost:9000"
    finally:
        stop_minio()
