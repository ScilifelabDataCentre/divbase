"""
e2e tests for the "divbase-cli query" commands

All tests are run against a docker compose setup with the entire DivBase stack running locally.

A project (CONSTANTS["QUERY_PROJECT"]) is made available with input files for the tests.
"""

import datetime
import logging
import os
import subprocess
import time
from contextlib import contextmanager
from pathlib import Path
from unittest.mock import patch

import boto3
import pytest
from celery import current_app
from sqlalchemy import select
from typer.testing import CliRunner

from divbase_api.exceptions import VCFDimensionsEntryMissingError
from divbase_api.models.task_history import CeleryTaskMeta, TaskHistoryDB
from divbase_api.services.queries import BcftoolsQueryManager
from divbase_api.services.s3_client import create_s3_file_manager
from divbase_api.services.task_history import _deserialize_celery_task_metadata
from divbase_api.worker.tasks import bcftools_pipe_task
from divbase_api.worker.worker_db import SyncSessionLocal
from divbase_cli.divbase_cli import app
from divbase_lib.api_schemas.task_history import TaskHistoryResult
from divbase_lib.exceptions import ProjectNotInConfigError

runner = CliRunner()


@pytest.fixture(autouse=True, scope="function")
def auto_clean_dimensions_entries_for_all_projects(clean_all_projects_dimensions):
    """Enable auto-cleanup of dimensions entries for all tests in this test file."""
    yield


def wait_for_task_complete(task_id: str, max_retries: int = 30) -> TaskHistoryResult:
    """
    For a given task_id, check the status of the task via the PostgreSQL results backend until it is complete or times out.
    Need to join CeleryTaskMeta with TaskHistoryDB to get timestamps for the mandatory fields of the pydantic model.
    Deserialization is needed because of how the celery results backed stores fields.
    """
    while max_retries > 0:
        with SyncSessionLocal() as db:
            stmt = (
                select(
                    *CeleryTaskMeta.__table__.c,
                    TaskHistoryDB.created_at,
                    TaskHistoryDB.started_at,
                    TaskHistoryDB.completed_at,
                )
                .join(TaskHistoryDB, CeleryTaskMeta.task_id == TaskHistoryDB.task_id)
                .where(CeleryTaskMeta.task_id == task_id)
            )

            result = db.execute(stmt).first()

            if result is None:
                time.sleep(1)
                max_retries -= 1
                continue

            task_dict = dict(result._mapping)

            status = task_dict.get("status")
            if status in ["SUCCESS", "FAILURE"]:
                deserialized = _deserialize_celery_task_metadata(task_dict)
                return TaskHistoryResult(**deserialized)

        time.sleep(1)
        max_retries -= 1

    pytest.fail(f"Task {task_id} didn't complete within timeout period")


@pytest.fixture(autouse=True)
def reset_query_projects_bucket(CONSTANTS):
    """Delete the merged.vcf.gz file from the query projects bucket before each test."""
    s3_resource = boto3.resource(
        "s3",
        endpoint_url=CONSTANTS["MINIO_URL"],
        aws_access_key_id=CONSTANTS["BAD_ACCESS_KEY"],
        aws_secret_access_key=CONSTANTS["BAD_SECRET_KEY"],
    )

    # pylance does not understand boto3 resource returns types, hence ignore below
    project_name = CONSTANTS["QUERY_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    bucket = s3_resource.Bucket(bucket_name)  # type: ignore
    bucket.object_versions.filter(Prefix="merged.vcf.gz").delete()
    yield


def test_sample_metadata_query(
    CONSTANTS,
    logged_in_edit_user_with_existing_config,
    run_update_dimensions,
    db_session_sync,
    project_map,
):
    """Test running a sample metadata query using the CLI."""
    project_name = CONSTANTS["QUERY_PROJECT"]
    project_id = project_map[project_name]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    run_update_dimensions(bucket_name=bucket_name, project_id=project_id, project_name=project_name)

    query_string = "Area:West of Ireland,Northern Portugal;Sex:F"
    expected_sample_ids = ["5a_HOM-I13", "5a_HOM-I14", "5a_HOM-I20", "5a_HOM-I21", "5a_HOM-I7", "1b_HOM-G58"]
    expected_filenames = ["HOM_20ind_17SNPs_last_10_samples.vcf.gz", "HOM_20ind_17SNPs_first_10_samples.vcf.gz"]

    command = f"query tsv '{query_string}' --project {project_name}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0

    assert query_string in result.stdout
    for sample_id in expected_sample_ids:
        assert sample_id in result.stdout
    for filename in expected_filenames:
        assert filename in result.stdout


def test_bcftools_pipe_query(
    CONSTANTS,
    logged_in_edit_user_with_existing_config,
    run_update_dimensions,
    db_session_sync,
    project_map,
):
    """Test running a bcftools pipe query using the CLI."""
    project_name = CONSTANTS["QUERY_PROJECT"]
    project_id = project_map[project_name]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    run_update_dimensions(bucket_name=bucket_name, project_id=project_id, project_name=project_name)
    tsv_filter = "Area:West of Ireland,Northern Portugal;"
    arg_command = "view -s SAMPLES; view -r 21:15000000-25000000"

    command = f"query bcftools-pipe --tsv-filter '{tsv_filter}' --command '{arg_command}' --project {project_name} "
    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert "Job submitted" in result.stdout

    task_id = result.stdout.strip().split()[-1]
    wait_for_task_complete(task_id=task_id)

    command = f"files list --project {project_name} "
    result = runner.invoke(app, command)

    assert result.exit_code == 0
    assert any("merged" in line and ".vcf.gz" in line for line in result.stdout.splitlines()), (
        "No merged VCF file found in output"
    )


def test_bcftools_pipe_fails_on_project_not_in_config(CONSTANTS, logged_in_edit_user_with_existing_config):
    project_name = "non_existent_project"
    tsv_filter = "Area:West of Ireland,Northern Portugal;"
    arg_command = "view -s SAMPLES"

    command = f"query bcftools-pipe --tsv-filter '{tsv_filter}' --command '{arg_command}' --project {project_name} "
    result = runner.invoke(app, command)
    assert isinstance(result.exception, ProjectNotInConfigError)


@pytest.mark.parametrize(
    "project_name,tsv_filter,command,expected_error",
    [
        # Malformed tsv filter (missing colon)
        ("DEFAULT", "Area West of Ireland", "DEFAULT", "SidecarInvalidFilterError"),
        # bad command
        ("DEFAULT", "DEFAULT", "invalid-command", "Unsupported bcftools command"),
        # empty command string
        ("DEFAULT", "DEFAULT", "", "Empty"),
    ],
)
def test_bcftools_pipe_query_errors(
    run_update_dimensions,
    db_session_sync,
    project_map,
    project_name,
    tsv_filter,
    command,
    expected_error,
    CONSTANTS,
    logged_in_edit_user_with_existing_config,
):
    """Test that validation errors cause task FAILURE status."""
    if "DEFAULT" in project_name:
        project_name = CONSTANTS["QUERY_PROJECT"]
    if "DEFAULT" in tsv_filter:
        tsv_filter = "Area:West of Ireland,Northern Portugal;"
    if "DEFAULT" in command:
        command = "view -s SAMPLES"

    project_id = project_map[project_name]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    run_update_dimensions(bucket_name=bucket_name, project_id=project_id, project_name=project_name)

    command_str = f"query bcftools-pipe --tsv-filter '{tsv_filter}' --command '{command}' --project {project_name} "
    response = runner.invoke(app, command_str)

    assert response.exit_code == 0
    task_id = response.stdout.strip().split()[-1]
    result = wait_for_task_complete(task_id=task_id)

    assert result.status == "FAILURE", f"Expected FAILURE status but got {result.status}"

    assert isinstance(result.result, dict), "Expected error result to be a dict"

    # Celery stores exceptions with these fields
    exc_type = str(result.result.get("exc_type", ""))
    exc_message = str(result.result.get("exc_message", ""))

    full_error = f"{exc_type} {exc_message}"
    assert expected_error in full_error, (
        f"Expected '{expected_error}' in error, got:\n"
        f"  exc_type: {exc_type}\n"
        f"  exc_message: {exc_message}\n"
        f"  full result: {result.result}"
    )


def test_get_task_status_by_task_id(CONSTANTS, logged_in_edit_user_with_existing_config, db_session_sync):
    """
    Get the status of a task by its ID.
    Uses the PostgreSQL Celery results backend to get task info.
    """
    project_name = CONSTANTS["QUERY_PROJECT"]
    tsv_filter = "Area:West of Ireland,Northern Portugal;"
    arg_command = "view -s SAMPLES; view -r 21:15000000-25000000"

    command = f"query bcftools-pipe --tsv-filter '{tsv_filter}' --command '{arg_command}' --project {project_name} "
    first_task_result = runner.invoke(app, command)
    assert first_task_result.exit_code == 0
    first_task_id = first_task_result.stdout.strip().split()[-1]

    second_task_result = runner.invoke(app, command)
    assert second_task_result.exit_code == 0
    second_task_id = second_task_result.stdout.strip().split()[-1]

    # Query the PostgreSQL results backend directly
    with SyncSessionLocal() as db:
        for task_id in [first_task_id, second_task_id]:
            stmt = select(CeleryTaskMeta).where(CeleryTaskMeta.task_id == task_id)
            result = db.execute(stmt).scalar_one_or_none()

            assert result is not None, f"Task {task_id} not found in results backend"
            assert result.task_id == task_id
            assert result.status in ["PENDING", "STARTED", "SUCCESS", "FAILURE"]


@patch(
    "divbase_api.worker.tasks.create_s3_file_manager",
    side_effect=lambda url=None: create_s3_file_manager(url="http://localhost:9002"),
)
def test_query_exits_when_vcf_file_version_is_outdated(
    mock_create_s3_manager,
    CONSTANTS,
    logged_in_edit_user_with_existing_config,
    fixtures_dir,
    run_update_dimensions,
    project_map,
    db_session_sync,
):
    """
    Test that updates the dimensions file, uploads a new version of a VCF file, then runs a query that should fail
    because the dimensions file expects an older version of the VCF file.
    """
    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    project_id = project_map[project_name]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    run_update_dimensions(bucket_name=bucket_name, project_id=project_id, project_name=project_name)

    def ensure_fixture_path(filename, fixture_dir="tests/fixtures"):
        if filename.startswith(fixture_dir):
            return filename
        return f"{fixture_dir}/{filename}"

    def patched_download_sample_metadata(metadata_tsv_name, bucket_name, s3_file_manager):
        """
        Patches the path for the sidecar metadata file so that it can be read from fixtures and not be downloaded.
        """
        return Path(ensure_fixture_path(metadata_tsv_name, fixture_dir="tests/fixtures"))

    def patched_download_vcf_files(files_to_download, bucket_name, s3_file_manager):
        """
        Needs the path in the worker container so that it is compatible with the docker exec patch below for running bcftools jobs.
        """
        pass

    with (
        patch("divbase_api.worker.tasks._download_sample_metadata", new=patched_download_sample_metadata),
        patch("divbase_api.worker.tasks._download_vcf_files", new=patched_download_vcf_files),
    ):
        test_file = (fixtures_dir / "HOM_20ind_17SNPs.1.vcf.gz").resolve()

        command = f"files upload {test_file}  --project {project_name} --disable-safe-mode"
        result = runner.invoke(app, command)

        assert result.exit_code == 0
        assert f"{str(test_file)}" in result.stdout

        params = {
            "tsv_filter": "Area:West of Ireland;Sex:F",
            "command": "view -s SAMPLES; view -r 1,4,6,21,24",
            "metadata_tsv_name": "sample_metadata_HOM_chr_split_version.tsv",
            "bucket_name": bucket_name,
            "user_name": "test-user",
            "project_id": project_id,
            "project_name": project_name,
        }
        with pytest.raises(ValueError) as excinfo:
            bcftools_pipe_task(**params)
        assert (
            "The VCF dimensions file is not up to date with the VCF files in the project. Please run 'divbase-cli dimensions update --project <project_name>' and then submit the query again."
            in str(excinfo.value)
        )


@pytest.mark.integration
@pytest.mark.parametrize(
    "params,expect_success,ensure_dimensions_file,expected_logs,expected_error_msgs",
    [
        # Case 0: expected to be fail, vcf dimensions file is empty so the early check for the file in the task raises an error
        (
            {
                "tsv_filter": "Area:West of Ireland;Sex:F",
                "command": "view -s SAMPLES; view -r 1,4,6,21,24",
                "metadata_tsv_name": "sample_metadata_HOM_chr_split_version.tsv",
                "project_name": "split-scaffold-project",
                "user_name": "test-user",
            },
            False,
            False,
            ["Starting bcftools_pipe_task"],
            [
                "The VCF dimensions index in project 'split-scaffold-project' is missing or empty. Please ensure that there are VCF files in the project and run:'divbase-cli dimensions update --project <project_name>'"
            ],
        ),
        # Case 1: expected to be sucessful, should lead to concat
        (
            {
                "tsv_filter": "Area:West of Ireland;Sex:F",
                "command": "view -s SAMPLES; view -r 1,4,6,21,24",
                "metadata_tsv_name": "sample_metadata_HOM_chr_split_version.tsv",
                "project_name": "split-scaffold-project",
                "user_name": "test-user",
            },
            True,
            True,
            [
                "Starting bcftools_pipe_task",
                "No unsupported sample sets found. Proceeding with bcftools pipeline.",
                "Sample names overlap between some temp files, will concat overlapping sets, then merge if needed and possible.",
                "Only one file remained after concatenation, renamed this file to",
                "Sorting the results file to ensure proper order of variants. Final results are in 'merged_",
                "bcftools processing completed successfully",
                "Cleaning up 14 temporary files",
            ],
            [],
        ),
        # Case 2: expected to fail since there are no scaffolds in the vcf files that match the -r query
        (
            {
                "tsv_filter": "Area:West of Ireland;Sex:F",
                "command": "view -s SAMPLES; view -r 31,34,36,321,324",
                "metadata_tsv_name": "sample_metadata_HOM_chr_split_version.tsv",
                "project_name": "split-scaffold-project",
                "user_name": "test-user",
                # project_id is added dynamically in the tests
            },
            False,
            True,
            [
                "Starting bcftools_pipe_task",
            ],
            [
                "Based on the 'view -r' query and the VCF scaffolds indexed in DivBase, there are no VCF files in the project that fulfill the query."
            ],
        ),
        # Case 3: expected to be sucessful, code should handle no tsv-filter in query
        (
            {
                "tsv_filter": "",
                "command": "view -s SAMPLES; view -r 1,4,6,21,24",
                "metadata_tsv_name": "sample_metadata_HOM_chr_split_version.tsv",
                "project_name": "split-scaffold-project",
                "user_name": "test-user",
                # project_id is added dynamically in the tests
            },
            True,
            True,
            [
                "Empty filter provided - returning ALL records. This may be a large result set.",
                "Starting bcftools_pipe_task",
                "'view -r' query requires scaffold '1'. It is present in file 'HOM_20ind_17SNPs.1.vcf.gz'",
                "'view -r' query requires scaffold '4'. It is present in file 'HOM_20ind_17SNPs.4.vcf.gz'",
                "'view -r' query requires scaffold '6'. It is present in file 'HOM_20ind_17SNPs.6.vcf.gz'",
                "'view -r' query requires scaffold '21'. It is present in file 'HOM_20ind_17SNPs.21.vcf.gz'",
                "'view -r' query requires scaffold '24'. It is present in file 'HOM_20ind_17SNPs.24.vcf.gz'",
                "No unsupported sample sets found. Proceeding with bcftools pipeline.",
                "Sample names overlap between some temp files, will concat overlapping sets, then merge if needed and possible.",
                "Only one file remained after concatenation, renamed this file to",
                "Sorting the results file to ensure proper order of variants. Final results are in 'merged_",
                "bcftools processing completed successfully",
            ],
            [],
        ),
        # Case 4: expected to be sucessful, should lead to merge
        (
            {
                "tsv_filter": "Area:Northern Portugal",
                "command": "view -s SAMPLES; view -r 21:15000000-25000000",
                "metadata_tsv_name": "sample_metadata.tsv",
                "project_name": "query-project",
                "user_name": "test-user",
                # project_id is added dynamically in the tests
            },
            True,
            True,
            [
                "Starting bcftools_pipe_task",
                "'view -r' query requires scaffold '21'. It is present in file 'HOM_20ind_17SNPs_last_10_samples.vcf.gz'.",
                "'view -r' query requires scaffold '21'. It is present in file 'HOM_20ind_17SNPs_first_10_samples.vcf.gz'.",
                "No unsupported sample sets found. Proceeding with bcftools pipeline.",
                "Sample names do not overlap between temp files, will continue with 'bcftools merge'",
                "Merged all temporary files into 'merged_unsorted_",
                "Sorting the results file to ensure proper order of variants. Final results are in 'merged_",
                "bcftools processing completed successfully",
                "Cleaning up 7 temporary files",
            ],
            [],
        ),
        # Case 5: expected to be sucessful, should lead one file being subset and renamed rather than bcftools merge/concat
        (
            {
                "tsv_filter": "Area:Northern Spanish shelf",
                "command": "view -s SAMPLES; view -r 1,4,6,21,24",
                "metadata_tsv_name": "sample_metadata_HOM_files_that_need_mixed_bcftools_concat_and_merge.tsv",
                "project_name": "mixed-concat-merge-project",
                "user_name": "test-user",
                # project_id is added dynamically in the tests
            },
            True,
            True,
            [
                "Starting bcftools_pipe_task",
                "'view -r' query requires scaffold '1'. It is present in file 'HOM_20ind_17SNPs.1.vcf.gz'.",
                "'view -r' query requires scaffold '4'. It is present in file 'HOM_20ind_17SNPs.4.vcf.gz'.",
                "'view -r' query requires scaffold '21'. It is present in file 'HOM_20ind_17SNPs.21.vcf.gz'.",
                "'view -r' query requires scaffold '1'. It is present in file 'HOM_20ind_17SNPs_changed_sample_names.vcf.gz'.",
                "'view -r' query requires scaffold '4'. It is present in file 'HOM_20ind_17SNPs_changed_sample_names.vcf.gz'.",
                "'view -r' query requires scaffold '6'. It is present in file 'HOM_20ind_17SNPs_changed_sample_names.vcf.gz'.",
                "'view -r' query requires scaffold '21'. It is present in file 'HOM_20ind_17SNPs_changed_sample_names.vcf.gz'.",
                "'view -r' query requires scaffold '24'. It is present in file 'HOM_20ind_17SNPs_changed_sample_names.vcf.gz'.",
                "No unsupported sample sets found. Proceeding with bcftools pipeline.",
                "Sample names overlap between some temp files, will concat overlapping sets, then merge if needed and possible.",
                "Merged all files (including concatenated files) into 'merged_unsorted_",
                "bcftools processing completed successfully",
                "Cleaning up 12 temporary files",
            ],
            [],
        ),
        # Case 6: expected to be sucessful, should lead to two different sample sets being concatenated, then merged with a third file
        (
            {
                "tsv_filter": "Area:Northern Spanish shelf,Iceland",
                "command": "view -s SAMPLES; view -r 1,4,6,8,13,18,21,24",
                "metadata_tsv_name": "sample_metadata_HOM_files_that_need_mixed_bcftools_concat_and_merge.tsv",
                "project_name": "mixed-concat-merge-project",
                "user_name": "test-user",
                # project_id is added dynamically in the tests
            },
            True,
            True,
            [
                "Starting bcftools_pipe_task",
                "'view -r' query requires scaffold '1'. It is present in file 'HOM_20ind_17SNPs.1.vcf.gz'.",
                "'view -r' query requires scaffold '4'. It is present in file 'HOM_20ind_17SNPs.4.vcf.gz'.",
                "'view -r' query requires scaffold '21'. It is present in file 'HOM_20ind_17SNPs.21.vcf.gz'.",
                "'view -r' query requires scaffold '1'. It is present in file 'HOM_20ind_17SNPs_changed_sample_names.vcf.gz'.",
                "'view -r' query requires scaffold '4'. It is present in file 'HOM_20ind_17SNPs_changed_sample_names.vcf.gz'.",
                "'view -r' query requires scaffold '6'. It is present in file 'HOM_20ind_17SNPs_changed_sample_names.vcf.gz'.",
                "'view -r' query requires scaffold '21'. It is present in file 'HOM_20ind_17SNPs_changed_sample_names.vcf.gz'.",
                "'view -r' query requires scaffold '24'. It is present in file 'HOM_20ind_17SNPs_changed_sample_names.vcf.gz'.",
                "'view -r' query requires scaffold '8'. It is present in file 'HOM_20ind_17SNPs.8_edit_new_sample_names.vcf.gz'.",
                "'view -r' query requires scaffold '13'. It is present in file 'HOM_20ind_17SNPs.13_edit_new_sample_names.vcf.gz'.",
                "'view -r' query requires scaffold '18'. It is present in file 'HOM_20ind_17SNPs.18_edit_new_sample_names.vcf.gz'.",
                "No unsupported sample sets found. Proceeding with bcftools pipeline.",
                "Sample names overlap between some temp files, will concat overlapping sets, then merge if needed and possible.",
                "Merged all files (including concatenated files) into 'merged_unsorted_",
                "bcftools processing completed successfully",
                "Cleaning up 19 temporary files",
            ],
            [],
        ),
    ],
)
def test_bcftools_pipe_cli_integration_with_eager_mode(
    tmp_path,
    CONSTANTS,
    caplog,
    params,
    expect_success,
    ensure_dimensions_file,
    expected_logs,
    expected_error_msgs,
    run_update_dimensions,
    project_map,
    db_session_sync,
):
    """
    This is a special integration test that allows for running bcftools-pipe queries directly in eager mode
    in a way that allows for catching the logs that otherwise would be printed inside the workers. For
    comparison, running a CLIrunner test will only give the "task submitted" log back.

    For this to work, a substantial amount of patching is needed. In short, since the task is run eagerly
    and directly, bcftools will need to be run with docker exec instead of subprocess (see BcftoolsQueryManager.run_bcftools).
    This is complicated by the fact that in the e2e process, files are transferred from the bucket to the
    worker. For the testing compose stack, the test buckets are built from ./tests/fixtures, so the workaround
    here is to patch out the transfer and have bcftools read all files directly from ./tests/fixtures. This
    works since the compose stack mounts ./tests/fixtures. However, for this test, the python code typically
    needs to look at ./tests/fixtures, but the docker exec needs to look at the mount in the container at
    /app/tests/fixtures. A lot of patching back and forth to ensure that every function looks in the right
    dir (locally or container) is thus needed.

    The benefit of all this patching is that now it is possible to parameterize the test for expected (worker) log outcomes!

    """
    project_name = params["project_name"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    params["bucket_name"] = bucket_name
    project_id = project_map[project_name]
    params["project_id"] = project_id

    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    caplog.set_level(logging.INFO)

    original_task_always_eager = current_app.conf.task_always_eager
    original_task_eager_propagates = current_app.conf.task_eager_propagates
    original_merge_or_concat_bcftools_temp_files = BcftoolsQueryManager.merge_or_concat_bcftools_temp_files
    original_boto3_client = boto3.client

    def ensure_fixture_path(filename, fixture_dir="tests/fixtures"):
        if filename.startswith(fixture_dir):
            return filename
        return f"{fixture_dir}/{filename}"

    def strip_fixture_dir(filename):
        fixture_dir = "tests/fixtures/"
        if filename.startswith(fixture_dir):
            return filename[len(fixture_dir) :]
        return filename

    def patched_download_sample_metadata(metadata_tsv_name, bucket_name, s3_file_manager):
        """
        Patches the path for the sidecar metadata file so that it can be read from fixtures and not be downloaded.
        """
        return Path(ensure_fixture_path(metadata_tsv_name, fixture_dir="tests/fixtures"))

    def patched_download_vcf_files(files_to_download, bucket_name, s3_file_manager):
        """
        Needs the path in the worker container so that it is compatible with the docker exec patch below for running bcftools jobs.
        """
        return [ensure_fixture_path(file_name, fixture_dir="/app/tests/fixtures") for file_name in files_to_download]

    def patched_run_bcftools(self, command: str) -> None:
        """
        Patches the working dir used when running bcftools commands inside the Docker container.
        """
        container_id = self.get_container_id(self.CONTAINER_NAME)
        logger = logging.getLogger("divbase_lib.queries")
        logger.debug(f"Executing command in container with ID: {container_id}")
        docker_cmd = ["docker", "exec", "-w", "/app/tests/fixtures", container_id, "bcftools"] + command.split()
        subprocess.run(docker_cmd, check=True)

    @contextmanager
    def patched_temp_file_management(self):
        """Context manager to handle temporary file cleanup, ensuring all temp files are in ./tests/fixtures."""
        self.temp_files = []
        try:
            yield self
        finally:
            if self.temp_files:
                logger = logging.getLogger("divbase_lib.queries")
                logger.info(f"Cleaning up {len(self.temp_files)} temporary files")
                temp_files_with_path = [ensure_fixture_path(f) for f in self.temp_files]
                self.cleanup_temp_files(temp_files_with_path)

    def patched_merge_or_concat_bcftools_temp_files(self, output_temp_files, identifier):
        """
        Patches a method that needs quite a bit of patching of submethods, hence nested patches...
        It ensures that all paths are correctly set to either ./tests/fixtures or /app/tests/fixtures depending
        on if python code or docker exec code needs to access the files.
        """
        original_rename = os.rename
        original_get_all_sample_names_from_vcf_files = self._get_all_sample_names_from_vcf_files
        original_group_vcfs_by_sample_set = self._group_vcfs_by_sample_set

        def patched_rename(src, dst):
            src = ensure_fixture_path(src)
            dst = ensure_fixture_path(dst)
            return original_rename(src, dst)

        def patched_get_all_sample_names_from_vcf_files(output_temp_files):
            output_temp_files = [ensure_fixture_path(f) for f in output_temp_files]
            return original_get_all_sample_names_from_vcf_files(output_temp_files)

        def patched_group_vcfs_by_sample_set(sample_names_per_VCF):
            stripped = {strip_fixture_dir(k): v for k, v in sample_names_per_VCF.items()}
            return original_group_vcfs_by_sample_set(stripped)

        with (
            patch("os.rename", new=patched_rename),
            patch.object(self, "_get_all_sample_names_from_vcf_files", new=patched_get_all_sample_names_from_vcf_files),
            patch.object(self, "_group_vcfs_by_sample_set", new=patched_group_vcfs_by_sample_set),
        ):
            return original_merge_or_concat_bcftools_temp_files(self, output_temp_files, identifier)

    def patched_upload_results_file(output_file, bucket_name, s3_file_manager):
        """
        Use the bucket_name from the test parameterization for uploading the results file
        """
        output_file = Path(ensure_fixture_path(str(output_file)))

        return s3_file_manager.upload_files(
            to_upload={output_file.name: output_file},
            bucket_name=bucket_name,
        )

    def patched_delete_job_files_from_worker(vcf_paths=None, metadata_path=None, output_file=None):
        """
        Only delete the output file, using the correct path. Don't delete the fixtures, since they should persist.
        """

        logger = logging.getLogger("divbase_api.worker.tasks")

        if output_file is not None:
            output_file = ensure_fixture_path(str(output_file))
            try:
                os.remove(output_file)
                logger.info(f"deleted {output_file}")
            except Exception as e:
                logger.warning(f"Could not delete output file from worker {output_file}: {e}")

    def patched_boto3_client(service_name, **kwargs):
        if service_name == "s3":
            kwargs["endpoint_url"] = CONSTANTS["MINIO_URL"]
        return original_boto3_client(service_name, **kwargs)

    def patched_prepare_txt_with_divbase_header_for_vcf(self, header_filename: str) -> None:
        """Create header file in the fixtures directory where the testing stack workers can find it"""
        header_path = ensure_fixture_path(header_filename)
        with open(header_path, "w") as file:
            file.write('##DivBase_created="This is a results file created by a DivBase query; ')
            file.write(f'Date={datetime.datetime.now().strftime("%a %b %d %H:%M:%S %Y")}"\n')

    if ensure_dimensions_file:
        run_update_dimensions(bucket_name=bucket_name, project_id=project_id, project_name=project_name)
    try:
        current_app.conf.update(
            task_always_eager=True,
            task_eager_propagates=True,
        )

        with (
            patch("boto3.client", side_effect=patched_boto3_client),
            patch("divbase_api.worker.tasks._download_sample_metadata", new=patched_download_sample_metadata),
            patch("divbase_api.services.queries.BcftoolsQueryManager.CONTAINER_NAME", "divbase-tests-worker-quick-1"),
            patch("divbase_api.worker.tasks._download_vcf_files", side_effect=patched_download_vcf_files),
            patch("divbase_api.services.queries.BcftoolsQueryManager.run_bcftools", new=patched_run_bcftools),
            patch(
                "divbase_api.services.queries.BcftoolsQueryManager.temp_file_management",
                new=patched_temp_file_management,
            ),
            patch(
                "divbase_api.services.queries.BcftoolsQueryManager.merge_or_concat_bcftools_temp_files",
                new=patched_merge_or_concat_bcftools_temp_files,
            ),
            patch("divbase_api.worker.tasks._upload_results_file", new=patched_upload_results_file),
            patch("divbase_api.worker.tasks._delete_job_files_from_worker", new=patched_delete_job_files_from_worker),
            patch(
                "divbase_api.services.queries.BcftoolsQueryManager._prepare_txt_with_divbase_header_for_vcf",
                new=patched_prepare_txt_with_divbase_header_for_vcf,
            ),
        ):
            if not expect_success:
                with pytest.raises((VCFDimensionsEntryMissingError, ValueError)) as e:
                    bcftools_pipe_task(**params)
                for msg in expected_error_msgs:
                    assert msg.replace("\n", "") in str(e.value).replace("\n", "")
            else:
                result = bcftools_pipe_task(**params)
                assert result is not None

            for log_msg in expected_logs:
                assert log_msg in caplog.text

            print(f"Captured logs:\n{caplog.text}")
    finally:
        current_app.conf.task_always_eager = original_task_always_eager
        current_app.conf.task_eager_propagates = original_task_eager_propagates
