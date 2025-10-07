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

from divbase_lib.s3_client import create_s3_file_manager
from divbase_worker.tasks import update_vcf_dimensions_task
from tests.helpers.docker_testing_stack_setup import start_compose_stack, stop_compose_stack
from tests.helpers.minio_setup import (
    MINIO_FAKE_ACCESS_KEY,
    MINIO_FAKE_SECRET_KEY,
    MINIO_URL,
    PROJECTS,
    setup_minio_data,
)

runner = CliRunner()

logger = logging.getLogger(__name__)


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
    }


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


@pytest.fixture
def run_update_dimensions(CONSTANTS):
    """
    Factory fixture that directly calls the update_vcf_dimensions_task task to create and update dimensions for the split-scaffold-project.
    Patches the S3 URL in the task to use the test MinIO url. Run directly to not have to poll for celery task completion.

    If run with test_minio_url, bucket_name = run_update_dimensions(), it uses the default bucket split-scaffold-project.
    It can also be run for other bucket names with: test_minio_url, bucket_name = run_update_dimensions("another-bucket-name")

    Since this fixture is not run in a celery worker, patch the delete_job_files_from_worker function just pass and do nothing.
    This ensures that no files are deleted and that no logging messages about deletion are printed.
    """
    test_minio_url = CONSTANTS["MINIO_URL"]
    default_bucket_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]

    def _run(bucket_name=default_bucket_name):
        with (
            patch("divbase_lib.vcf_dimension_indexing.logger.info") as mock_info,
            patch("divbase_worker.tasks.create_s3_file_manager") as mock_create_s3_manager,
            patch("divbase_worker.tasks.delete_job_files_from_worker") as mock_delete_job_files,
        ):

            def append_test_fixture_info_to_log(msg, *args, **kwargs):
                if "No VCF dimensions file found in the bucket:" in msg:
                    msg = "LOG CALL BY TEST FIXTURE (run_update_dimensions): " + msg
                return mock_info.original(msg, *args, **kwargs)

            def patched_delete_job_files_from_worker(vcf_paths=None, metadata_path=None, output_file=None):
                pass

            mock_info.original = logger.info
            mock_info.side_effect = append_test_fixture_info_to_log
            mock_create_s3_manager.side_effect = lambda url=None: create_s3_file_manager(url=test_minio_url)
            mock_delete_job_files.side_effect = patched_delete_job_files_from_worker
            result = update_vcf_dimensions_task(bucket_name=bucket_name)
            assert result["status"] == "completed"

            # Workaround for garbage collection since patching in tmp_path across all layers turned out to be very complex...
            files_indexed_by_this_job = result.get("VCF files that were added to dimensions file by this job", [])
            for file_path in files_indexed_by_this_job:
                file_path_obj = Path(file_path).resolve()
                fixtures_dir = (Path(__file__).parent.parent / "fixtures").resolve()
                if fixtures_dir not in file_path_obj.parents:
                    file_path_obj.unlink()

    return _run
