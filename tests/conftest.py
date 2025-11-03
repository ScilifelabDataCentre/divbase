"""
Top-level pytest configuration for the DivBase project.
It handles spinning up the job system docker stack for the duration of the test session, and the tear-down afterwards.

It also collects fixtures and constants that are needed across multiple test modules. Storing them here in the top-level ensures
that the imports work correctly.
"""

import logging
import os
from pathlib import Path
from unittest.mock import patch

import httpx
import pytest
from typer.testing import CliRunner

from divbase_api.worker.tasks import update_vcf_dimensions_task
from divbase_lib.s3_client import create_s3_file_manager
from divbase_lib.vcf_dimension_indexing import create_vcf_dimension_manager
from tests.helpers.api_setup import ADMIN_CREDENTIALS, TEST_USERS, setup_api_data
from tests.helpers.docker_testing_stack_setup import start_compose_stack, stop_compose_stack
from tests.helpers.minio_setup import (
    MINIO_FAKE_ACCESS_KEY,
    MINIO_FAKE_SECRET_KEY,
    MINIO_URL,
    PROJECTS,
    setup_minio_data,
)

# Set env vars before the fixtures to ensure that they are available at the right time
os.environ["CELERY_BROKER_URL"] = "pyamqp://guest@localhost:5673//"
os.environ["CELERY_RESULT_BACKEND"] = "redis://localhost:6380/0"
os.environ["FLOWER_USER"] = "floweradmin"
os.environ["FLOWER_PASSWORD"] = "badpassword"
os.environ["FLOWER_BASE_URL"] = "http://localhost:5556"

os.environ["DIVBASE_S3_ACCESS_KEY"] = MINIO_FAKE_ACCESS_KEY
os.environ["DIVBASE_S3_SECRET_KEY"] = MINIO_FAKE_SECRET_KEY

os.environ["DATABASE_URL"] = "postgresql+asyncpg://divbase_user:badpassword@postgres:5432/divbase_db"
os.environ["WORKER_SERVICE_EMAIL"] = "worker@divbase.com"
os.environ["WORKER_SERVICE_PASSWORD"] = "badpassword"


api_base_url = os.environ["DIVBASE_API_URL"]

runner = CliRunner()

logger = logging.getLogger(__name__)


@pytest.fixture(autouse=True, scope="function")
def clean_tmp_config_token_dir():
    """
    To avoid test pollution, ensure that any config or token files created in tests
    are removed after each test function and are not created where the local dev/user has their config or token file.

    Related to this fixture is cli_config.py where the cli_settings instance is created at module load time.
    Using env variables we point to these "testing" paths instead for the config and token paths.
    """
    test_config_path = Path("tests/fixtures/tmp/config.yaml")
    test_token_path = Path("tests/fixtures/tmp/.fakesecrets")

    test_config_path.unlink(missing_ok=True)
    test_token_path.unlink(missing_ok=True)
    yield
    test_config_path.unlink(missing_ok=True)
    test_token_path.unlink(missing_ok=True)


@pytest.fixture(scope="session")
def CONSTANTS():
    return {
        "BAD_ACCESS_KEY": MINIO_FAKE_ACCESS_KEY,
        "BAD_SECRET_KEY": MINIO_FAKE_SECRET_KEY,
        "MINIO_URL": MINIO_URL,
        "DEFAULT_PROJECT": "project1",
        "NON_DEFAULT_PROJECT": "project2",
        "QUERY_PROJECT": "query-project",
        "SPLIT_SCAFFOLD_PROJECT": "split-scaffold-project",
        "CLEANED_PROJECT": "cleaned-project",
        "EMPTY_PROJECT": "empty-project",
        "PROJECT_CONTENTS": PROJECTS,
        "FILES_TO_UPLOAD_DOWNLOAD": ["file1.txt", "file2.txt", "file3.txt"],
        "ADMIN_CREDENTIALS": ADMIN_CREDENTIALS,
        "TEST_USERS": TEST_USERS,
    }


@pytest.fixture(autouse=True, scope="session")
def docker_testing_stack():
    """
    Start job system docker stack, and stop after all tests run.
    """
    try:
        start_compose_stack()
        setup_minio_data()
        setup_api_data()
        yield
    finally:
        stop_compose_stack()


@pytest.fixture
def run_update_dimensions(CONSTANTS):
    """
    Factory fixture that directly calls the update_vcf_dimensions_task task to create and update dimensions for the split-scaffold-project.
    Uses VCFDimensionIndexManager for cleanup/setup.
    """
    default_bucket_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]

    def _run(bucket_name=default_bucket_name):
        login_response = httpx.post(
            f"{api_base_url}/v1/auth/login",
            data={
                "grant_type": "password",
                "username": os.environ["WORKER_SERVICE_EMAIL"],
                "password": os.environ["WORKER_SERVICE_PASSWORD"],
            },
            headers={"Content-Type": "application/x-www-form-urlencoded"},
        )
        login_response.raise_for_status()
        auth_token = login_response.json()["access_token"]

        # Step 2: Use VCFDimensionIndexManager to clean up all VCF metadata and skipped files
        dimension_manager = create_vcf_dimension_manager(bucket_name=bucket_name, auth_token=auth_token)
        try:
            vcf_dimensions_data = dimension_manager.get_dimensions_info()
            project_id = dimension_manager._project_id

            for entry in vcf_dimensions_data.get("vcf_files", []):
                vcf_file = entry["vcf_file_s3_key"]
                try:
                    dimension_manager.delete_vcf_metadata(vcf_file, project_id)
                except Exception as e:
                    print(f"Warning: Failed to delete VCF metadata for {vcf_file}: {e}")

            skipped_files = dimension_manager.get_skipped_files()
            for vcf_file in skipped_files:
                try:
                    dimension_manager.delete_skipped_vcf(vcf_file, project_id)
                except Exception as e:
                    print(f"Warning: Failed to delete skipped VCF entry for {vcf_file}: {e}")

        except Exception as e:
            print(f"Error cleaning up dimensions for {bucket_name}: {e}")

        with patch("divbase_api.worker.tasks.create_s3_file_manager") as mock_create_s3_manager:
            mock_create_s3_manager.side_effect = lambda url=None: create_s3_file_manager(url=CONSTANTS["MINIO_URL"])
            result = update_vcf_dimensions_task(bucket_name=bucket_name)
        return result

    return _run
