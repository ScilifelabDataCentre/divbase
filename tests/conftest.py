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

import pytest
from typer.testing import CliRunner

from divbase_api.worker.crud_dimensions import (
    delete_skipped_vcf,
    delete_vcf_metadata,
    get_skipped_vcfs_by_project_worker,
    get_vcf_metadata_by_project,
)
from divbase_api.worker.tasks import update_vcf_dimensions_task
from divbase_api.worker.worker_db import SyncSessionLocal
from divbase_lib.s3_client import create_s3_file_manager
from tests.helpers.api_setup import (
    ADMIN_CREDENTIALS,
    TEST_USERS,
    get_project_map,
    setup_api_data,
)
from tests.helpers.docker_testing_stack_setup import start_compose_stack, stop_compose_stack
from tests.helpers.minio_setup import (
    MINIO_FAKE_ACCESS_KEY,
    MINIO_FAKE_SECRET_KEY,
    MINIO_URL,
    PROJECTS,
    setup_minio_data,
)

os.environ["DIVBASE_S3_ACCESS_KEY"] = MINIO_FAKE_ACCESS_KEY
os.environ["DIVBASE_S3_SECRET_KEY"] = MINIO_FAKE_SECRET_KEY

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


@pytest.fixture(scope="session")
def project_map(docker_testing_stack):
    """Get project_map after API setup."""
    return get_project_map()


@pytest.fixture
def db_session_sync():
    with SyncSessionLocal() as db:
        yield db


@pytest.fixture
def clean_vcf_dimensions():
    """
    Factory fixture to clean VCF dimensions for a specific project.
    Usage: clean_vcf_dimensions(db_session_sync, project_id)
    """

    # TODO perhaps there could be a a delete all per project function in the crud_dimensions module?
    def _clean(db, project_id):
        try:
            vcf_dimensions_data = get_vcf_metadata_by_project(project_id=project_id, db=db)

            for entry in vcf_dimensions_data.get("vcf_files", []):
                vcf_file = entry["vcf_file_s3_key"]
                try:
                    delete_vcf_metadata(db=db, vcf_file_s3_key=vcf_file, project_id=project_id)
                except Exception as e:
                    print(f"Warning: Failed to delete VCF metadata for {vcf_file}: {e}")

            skipped_files = get_skipped_vcfs_by_project_worker(db=db, project_id=project_id)
            for vcf_file in skipped_files:
                try:
                    delete_skipped_vcf(db=db, vcf_file_s3_key=vcf_file, project_id=project_id)
                except Exception as e:
                    print(f"Warning: Failed to delete skipped VCF entry for {vcf_file}: {e}")

        except Exception as e:
            print(f"Error cleaning up dimensions for project {project_id}: {e}")

    return _clean


@pytest.fixture
def run_update_dimensions(CONSTANTS):
    """
    Factory fixture that runs update_vcf_dimensions_task.
    Usage: run_update_dimensions(bucket_name)
    """

    def _update(bucket_name="split-scaffold-project", project_id=None, user_name="Test User"):
        with patch("divbase_api.worker.tasks.create_s3_file_manager") as mock_create_s3_manager:
            mock_create_s3_manager.side_effect = lambda url=None: create_s3_file_manager(url=CONSTANTS["MINIO_URL"])
            result = update_vcf_dimensions_task(bucket_name=bucket_name, project_id=project_id, user_name=user_name)
        return result

    return _update
