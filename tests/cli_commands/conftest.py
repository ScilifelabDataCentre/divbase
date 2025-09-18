"""
Setup for pytest fixtures for e2e testing of the CLI commands.

Pytest fixtures set up (and tear down) a test environment with the full DivBase stack running locally,

The S3FileManager class is patched in all tests to use this test MinIO server,
it is autoused, so it does not need to be specified in each test.
"""

import logging
from pathlib import Path
from unittest.mock import patch

import pytest
from typer.testing import CliRunner

from divbase_cli.divbase_cli import app
from divbase_lib.s3_client import create_s3_file_manager
from divbase_worker.tasks import update_vcf_dimensions_task
from tests.helpers.minio_setup import (
    MINIO_FAKE_ACCESS_KEY,
    MINIO_FAKE_SECRET_KEY,
    MINIO_URL,
    PROJECTS,
)

runner = CliRunner()

logger = logging.getLogger(__name__)


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


@pytest.fixture
def tmp_config_path(tmp_path):
    """
    Fixture to provide a path to where the a configuration can be created.
    """
    return tmp_path / "test_config.yaml"


@pytest.fixture
def fresh_config(tmp_path):
    """
    Fixture to provide a path to a pre-existing configuration file.
    """
    fresh_config_path = tmp_path / ".divbase_config.yaml"
    create_command = f"config create --config-file {fresh_config_path}"
    result = runner.invoke(app, create_command)

    assert result.exit_code == 0
    assert tmp_path.exists(), "Config file was not created"

    return fresh_config_path


@pytest.fixture
def user_config_path(tmp_path, CONSTANTS):
    """
    Fixture to provide a path to an "existing" user configuration file with
    some existing projects and a default project set.
    """
    existing_config_path = tmp_path / ".divbase_config.yaml"
    create_command = f"config create --config-file {existing_config_path}"
    result = runner.invoke(app, create_command)
    assert result.exit_code == 0

    for project in CONSTANTS["PROJECT_CONTENTS"]:
        add_command = f"config add-project {project} --divbase-url http://localhost:8001 --s3-url {MINIO_URL} --config {existing_config_path}"
        result = runner.invoke(app, add_command)
        assert result.exit_code == 0

    set_default_command = f"config set-default {CONSTANTS['DEFAULT_PROJECT']} --config {existing_config_path}"
    result = runner.invoke(app, set_default_command)
    assert result.exit_code == 0

    return existing_config_path


@pytest.fixture
def fixtures_dir():
    """Path to the fixtures directory."""
    return Path(__file__).parent.parent / "fixtures"


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
