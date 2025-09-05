"""
e2e tests for the "divbase-cli query" commands

All tests are run against a docker compose setup with the entire DivBase stack running locally.

A project (CONSTANTS["QUERY_PROJECT"]) is made available with input files for the tests.
"""

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
from typer.testing import CliRunner

from divbase_tools.divbase_cli import app
from divbase_tools.exceptions import ProjectNotInConfigError
from divbase_tools.queries import BcftoolsQueryManager
from divbase_tools.tasks import bcftools_pipe_task
from tests.helpers.minio_setup import MINIO_URL

runner = CliRunner()


def wait_for_task_complete(task_id: str, config_file: str, max_retries: int = 30) -> None:
    """Given a task_id, check the status of the task via the CLI until it is complete or times out."""
    command = f"query task-status {task_id} --config {config_file}"
    while max_retries > 0:
        result = runner.invoke(app, command)
        if (
            "FAILURE" in result.stdout
            or "SUCCESS" in result.stdout
            or "'status': 'completed'" in result.stdout
            or "completed" in result.stdout
            or "FAIL" in result.stdout
        ):
            return
        time.sleep(1)
        max_retries -= 1
    pytest.fail(f"Task didn't complete retries. Last status: {result.stdout}")


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
    bucket = s3_resource.Bucket(CONSTANTS["QUERY_PROJECT"])  # type: ignore
    bucket.object_versions.filter(Prefix="merged.vcf.gz").delete()
    yield


def test_sample_metadata_query(CONSTANTS, user_config_path):
    """Test running a sample metadata query using the CLI."""
    project_name = CONSTANTS["QUERY_PROJECT"]
    query_string = "Area:West of Ireland,Northern Portugal;Sex:F"
    expected_sample_ids = "['5a_HOM-I13', '5a_HOM-I14', '5a_HOM-I20', '5a_HOM-I21', '5a_HOM-I7', '1b_HOM-G58']"
    expected_filenames = "['HOM_20ind_17SNPs_last_10_samples.vcf.gz', 'HOM_20ind_17SNPs_first_10_samples.vcf.gz']"

    command = f"query tsv '{query_string}' --project {project_name} --config {user_config_path}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0

    assert query_string in result.stdout
    assert expected_sample_ids in result.stdout
    assert expected_filenames in result.stdout


def test_bcftools_pipe_query(user_config_path, CONSTANTS):
    """Test running a bcftools pipe query using the CLI."""
    project_name = CONSTANTS["QUERY_PROJECT"]
    tsv_filter = "Area:West of Ireland,Northern Portugal;"
    arg_command = "view -s SAMPLES; view -r 21:15000000-25000000"

    command = f"query bcftools-pipe --tsv-filter '{tsv_filter}' --command '{arg_command}' --project {project_name} --config {user_config_path}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert "Job submitted" in result.stdout

    task_id = result.stdout.strip().split()[-1]
    wait_for_task_complete(task_id=task_id, config_file=user_config_path)

    command = f"files list --project {project_name} --config {user_config_path}"
    result = runner.invoke(app, command)

    assert result.exit_code == 0
    assert any("merged" in line and ".vcf.gz" in line for line in result.stdout.splitlines()), (
        "No merged VCF file found in output"
    )


def test_bcftools_pipe_fails_on_project_not_in_config(CONSTANTS, user_config_path):
    project_name = "non_existent_project"
    tsv_filter = "Area:West of Ireland,Northern Portugal;"
    arg_command = "view -s SAMPLES"

    command = f"query bcftools-pipe --tsv-filter '{tsv_filter}' --command '{arg_command}' --project {project_name} --config {user_config_path}"
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
def test_bcftools_pipe_query_errors(project_name, tsv_filter, command, expected_error, CONSTANTS, user_config_path):
    """
    Test bad formatted input raises errors

    TODO - these sorts of errors in the future should be handled by the API, and not sent to Celery.
    This test will need to be rewritten then, hence why it is quite sloppy now.
    """
    if "DEFAULT" in project_name:
        project_name = CONSTANTS["QUERY_PROJECT"]
    if "DEFAULT" in tsv_filter:
        tsv_filter = "Area:West of Ireland,Northern Portugal;"
    if "DEFAULT" in command:
        command = "view -s SAMPLES"

    command = f"query bcftools-pipe --tsv-filter '{tsv_filter}' --command '{command}' --project {project_name} --config {user_config_path}"
    result = runner.invoke(app, command)

    task_id = result.stdout.strip().split()[-1]
    # TODO, assertion below should become 1, when API has the role of validating input.
    assert result.exit_code == 0
    wait_for_task_complete(task_id=task_id, config_file=user_config_path)

    command = f"query task-status {task_id} --config {user_config_path}"
    job_status = runner.invoke(app, command)
    assert job_status.exit_code == 0
    assert expected_error in job_status.stdout


def test_get_task_status_by_task_id(CONSTANTS, user_config_path):
    """Get the status of a task by its ID."""
    project_name = CONSTANTS["QUERY_PROJECT"]
    tsv_filter = "Area:West of Ireland,Northern Portugal;"
    arg_command = "view -s SAMPLES; view -r 21:15000000-25000000"

    command = f"query bcftools-pipe --tsv-filter '{tsv_filter}' --command '{arg_command}' --project {project_name} --config {user_config_path}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0
    task_id = result.stdout.strip().split()[-1]

    other_result = runner.invoke(app, command)
    assert other_result.exit_code == 0
    other_task_id = other_result.stdout.strip().split()[-1]

    command = f"query task-status {task_id} --config {user_config_path}"
    result = runner.invoke(app, command)
    print(result.stdout)
    assert result.exit_code == 0
    assert task_id in result.stdout
    assert other_task_id not in result.stdout


@pytest.mark.integration
@pytest.mark.parametrize(
    "params,expect_success,ensure_dimensions_file,expected_logs,expected_error_msgs",
    [
        # Case: expected to be fail, vcf dimensions file is empty
        (
            {
                "tsv_filter": "Area:West of Ireland;Sex:F",
                "command": "view -s SAMPLES; view -r 1,4,6,21,24",
                "metadata_tsv_name": "sample_metadata_HOM_chr_split_version.tsv",
                "bucket_name": "split-scaffold-project",
                "user_name": "test-user",
            },
            True,
            False,
            [
                "Starting bcftools_pipe_task",
                "Sample names overlap between some temp files, will concat overlapping sets, then merge if needed and possible.",
            ],
            [],
        ),
        # case: expected to be sucessful, should lead to concat
        (
            {
                "tsv_filter": "Area:West of Ireland;Sex:F",
                "command": "view -s SAMPLES; view -r 1,4,6,21,24",
                "metadata_tsv_name": "sample_metadata_HOM_chr_split_version.tsv",
                "bucket_name": "split-scaffold-project",
                "user_name": "test-user",
            },
            True,
            True,
            [
                "Starting bcftools_pipe_task",
                "Sample names overlap between some temp files, will concat overlapping sets, then merge if needed and possible.",
            ],
            [],
        ),
        # case
        (
            {
                "tsv_filter": "Area:West of Ireland;Sex:F",
                "command": "view -s SAMPLES; view -r 31,34,36,321,324",
                "metadata_tsv_name": "sample_metadata_HOM_chr_split_version.tsv",
                "bucket_name": "split-scaffold-project",
                "user_name": "test-user",
            },
            False,
            True,
            [
                "Starting bcftools_pipe_task",
            ],
            [
                "Based on the 'view -r' query and the VCF scaffolds indexed in DivBase, there are no VCF files in the project that fulfills the query. Please try another -r query with scaffolds/chromosomes that are present in the VCF files.To see a list of all unique scaffolds that are present across the VCF files in the project:'DIVBASE_ENV=local divbase-cli dimensions show --unique-scaffolds --project <PROJECT_NAME>"
            ],
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
    delete_dimensions_file_from_a_bucket,
):
    """
    This is a special integration test that allows for running bcftools-pipe queries
    directly in eager mode in a way that allows for catching the logs that otherwise
    would be printed inside the workers. For comparison, running a CLIrunner test will only
    give the "task submitted" log back.

    For this to work, a substantial amount of patching is needed. In short, since the task is
    run eagerly and directly, bcftools will need to be run with docker exec instead of subprocess
    (see BcftoolsQueryManager.run_bcftools). This is complicated by the fact that in the e2e process,
    files are transferred from the bucket to the worker. For the testing compose stack, the test buckets
    are built from ./tests/fixtures, so the workaround here is to patch out the transfer and have bcftools
    read all files directly from ./tests/fixtures. This works since the compose stack mounts ./tests/fixtures.
    However, for this test, the python code typically needs to look at ./tests/fixtures, but the docker exec
    needs to look at the mount in the container at /app/tests/fixtures. A lot of patching back and forth
    to ensure that every function looks in the right dir (locally or container) is thus needed.

    The benefit of all this patching is that now it is possible to parameterize the test for expected (worker) log outcomes!

    """
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
        logger = logging.getLogger("divbase_tools.queries")
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
                logger = logging.getLogger("divbase_tools.queries")
                logger.info(f"Cleaning up {len(self.temp_files)} temporary files")
                temp_files_with_path = [ensure_fixture_path(f) for f in self.temp_files]
                self.cleanup_temp_files(temp_files_with_path)

    def patched_merge_or_concat_bcftools_temp_files(self, output_temp_files, identifier):
        """
        Patches a method that needs quite a bit of patching of submethods, hence nested patches...
        """
        output_temp_files = [ensure_fixture_path(f) for f in output_temp_files]

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

    def patched_delete_job_files_from_worker(vcf_paths, metadata_path=None, output_file=None):
        """
        Only delete the output file, using the correct path. Don't delete the fixtures, since they should persist.
        """

        logger = logging.getLogger("divbase_tools.tasks")

        if output_file is not None:
            output_file = ensure_fixture_path(str(output_file))
            try:
                os.remove(output_file)
                logger.info(f"deleted {output_file}")
            except Exception as e:
                logger.warning(f"Could not delete output file from worker {output_file}: {e}")

    def patched_boto3_client(service_name, **kwargs):
        if service_name == "s3":
            kwargs["endpoint_url"] = MINIO_URL
        return original_boto3_client(service_name, **kwargs)

    if ensure_dimensions_file:
        run_update_dimensions(bucket_name=params["bucket_name"])
    else:
        delete_dimensions_file_from_a_bucket(params["bucket_name"])

    try:
        current_app.conf.update(
            task_always_eager=True,
            task_eager_propagates=True,
        )

        with (
            patch("boto3.client", side_effect=patched_boto3_client),
            patch("divbase_tools.tasks.download_sample_metadata", new=patched_download_sample_metadata),
            patch("divbase_tools.queries.BcftoolsQueryManager.CONTAINER_NAME", "divbase-tests-worker-quick-1"),
            patch("divbase_tools.tasks.download_vcf_files", side_effect=patched_download_vcf_files),
            patch("divbase_tools.queries.BcftoolsQueryManager.run_bcftools", new=patched_run_bcftools),
            patch("divbase_tools.queries.BcftoolsQueryManager.temp_file_management", new=patched_temp_file_management),
            patch(
                "divbase_tools.queries.BcftoolsQueryManager.merge_or_concat_bcftools_temp_files",
                new=patched_merge_or_concat_bcftools_temp_files,
            ),
            patch("divbase_tools.tasks.upload_results_file", new=patched_upload_results_file),
            patch("divbase_tools.tasks.delete_job_files_from_worker", new=patched_delete_job_files_from_worker),
        ):
            if not expect_success:
                with pytest.raises(ValueError) as excinfo:
                    bcftools_pipe_task(**params)
                for msg in expected_error_msgs:
                    assert msg.replace("\n", "") in str(excinfo.value).replace("\n", "")
            else:
                result = bcftools_pipe_task(**params)
                assert result is not None

            for log_msg in expected_logs:
                assert log_msg in caplog.text

            print(f"Captured logs:\n{caplog.text}")
    finally:
        current_app.conf.task_always_eager = original_task_always_eager
        current_app.conf.task_eager_propagates = original_task_eager_propagates
