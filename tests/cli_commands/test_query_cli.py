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
import httpx
import pytest
import yaml
from celery import current_app
from typer.testing import CliRunner

from divbase_cli.divbase_cli import app
from divbase_lib.exceptions import ProjectNotInConfigError
from divbase_lib.queries import BcftoolsQueryManager
from divbase_lib.s3_client import create_s3_file_manager
from divbase_lib.vcf_dimension_indexing import DIMENSIONS_FILE_NAME
from divbase_worker.tasks import bcftools_pipe_task
from tests.helpers.minio_setup import MINIO_URL

runner = CliRunner()


@pytest.fixture(autouse=True)
def clean_dimensions(user_config_path, CONSTANTS):
    """
    Remove the dimensions file and any merged output files before each test.
    Used in all tests in this module.
    """
    for project_name in CONSTANTS["PROJECT_CONTENTS"]:
        s3_client = boto3.client(
            "s3",
            endpoint_url=CONSTANTS["MINIO_URL"],
            aws_access_key_id=CONSTANTS["BAD_ACCESS_KEY"],
            aws_secret_access_key=CONSTANTS["BAD_SECRET_KEY"],
        )
        s3_client.delete_object(Bucket=project_name, Key=DIMENSIONS_FILE_NAME)

        try:
            # List all objects and delete any that match merged*.vcf.gz pattern, which caused contamination issues in parameterized tests)
            response = s3_client.list_objects_v2(Bucket=project_name)
            if "Contents" in response:
                for obj in response["Contents"]:
                    if obj["Key"].startswith("merged_") and obj["Key"].endswith(".vcf.gz"):
                        s3_client.delete_object(Bucket=project_name, Key=obj["Key"])
                        print(f"Deleted lingering merged file: {obj['Key']} from {project_name}")
        except Exception as e:
            print(f"Error cleaning up merged files: {e}")

        time.sleep(1)  # Give MinIO/S3 time to propagate deletion

    yield


def wait_for_task_complete(task_id: str, max_retries: int = 30):
    """Given a task_id, check the status of the task via the CLI until it is complete or times out."""
    command = f"query task-status {task_id}"
    while max_retries > 0:
        result = runner.invoke(app, command)
        # Add checks for error keywords
        if (
            "FAILURE" in result.stdout
            or "SUCCESS" in result.stdout
            or "'status': 'completed'" in result.stdout
            or "completed" in result.stdout
            or "FAIL" in result.stdout
            or "SidecarInvalidFilterError" in result.stdout
            or "Unsupported bcftools command" in result.stdout
            or "Empty command provided" in result.stdout
            or "'status': 'error'" in result.stdout
        ):
            return result
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


def test_bcftools_pipe_query(run_update_dimensions, user_config_path, CONSTANTS):
    """Test running a bcftools pipe query using the CLI."""
    project_name = CONSTANTS["QUERY_PROJECT"]
    tsv_filter = "Area:West of Ireland,Northern Portugal;"
    arg_command = "view -s SAMPLES; view -r 21:15000000-25000000"
    run_update_dimensions(bucket_name=project_name)

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


def test_bcftools_pipe_fails_on_project_not_in_config(CONSTANTS, user_config_path):
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
    run_update_dimensions, project_name, tsv_filter, command, expected_error, CONSTANTS, user_config_path
):
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
    run_update_dimensions(bucket_name=project_name)

    command = f"query bcftools-pipe --tsv-filter '{tsv_filter}' --command '{command}' --project {project_name} "
    result = runner.invoke(app, command)

    task_id = result.stdout.strip().split()[-1]
    # TODO, assertion below should become 1, when API has the role of validating input.
    assert result.exit_code == 0
    result = wait_for_task_complete(task_id=task_id)
    assert expected_error in result.stdout


def test_get_task_status_by_task_id(CONSTANTS, user_config_path):
    """Get the status of a task by its ID. Uses flower API via get_task_history to get the task info.
    Note that this does not test the CLI command for testing task status.
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

    flower_user = os.environ["FLOWER_USER"]
    flower_password = os.environ["FLOWER_PASSWORD"]
    flower_base_url = os.environ["FLOWER_BASE_URL"]

    for task_id in [first_task_id, second_task_id]:
        flower_url = f"{flower_base_url}/api/task/info/{task_id}"
        auth = (flower_user, flower_password)
        response = httpx.get(flower_url, auth=auth, timeout=3.0)
        tasks_status = response.json()
        assert tasks_status.get("uuid") == task_id


@pytest.mark.integration
@pytest.mark.parametrize(
    "params,expect_success,ensure_dimensions_file,expected_logs,expected_error_msgs",
    [
        # Case: expected to be fail, vcf dimensions file is empty so the early check for the file in the task raises an error
        (
            {
                "tsv_filter": "Area:West of Ireland;Sex:F",
                "command": "view -s SAMPLES; view -r 1,4,6,21,24",
                "metadata_tsv_name": "sample_metadata_HOM_chr_split_version.tsv",
                "bucket_name": "split-scaffold-project",
                "user_name": "test-user",
            },
            False,
            False,
            [
                "Starting bcftools_pipe_task",
                "No VCF dimensions file found in the bucket: split-scaffold-project.",
            ],
            ["The VCF dimensions file in project split-scaffold-project is missing or empty."],
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
                "No unsupported sample sets found. Proceeding with bcftools pipeline.",
                "Sample names overlap between some temp files, will concat overlapping sets, then merge if needed and possible.",
                "Only one file remained after concatenation, renamed this file to",
                "Sorting the results file to ensure proper order of variants. Final results are in 'merged_",
                "bcftools processing completed successfully",
                "Cleaning up 12 temporary files",
            ],
            [],
        ),
        # case: expected to fail since there are no scaffolds in the vcf files that match the -r query
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
        # case: expected to be sucessful, code should handle no tsv-filter in query
        (
            {
                "tsv_filter": "",
                "command": "view -s SAMPLES; view -r 1,4,6,21,24",
                "metadata_tsv_name": "sample_metadata_HOM_chr_split_version.tsv",
                "bucket_name": "split-scaffold-project",
                "user_name": "test-user",
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
        # case: expected to be sucessful, should lead to merge
        (
            {
                "tsv_filter": "Area:Northern Portugal",
                "command": "view -s SAMPLES; view -r 21:15000000-25000000",
                "metadata_tsv_name": "sample_metadata.tsv",
                "bucket_name": "query-project",
                "user_name": "test-user",
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
                "Cleaning up 5 temporary files",
            ],
            [],
        ),
        # case: expected to be sucessful, should lead one file being subset and renamed rather than bcftools merge/concat
        (
            {
                "tsv_filter": "Area:Northern Spanish shelf",
                "command": "view -s SAMPLES; view -r 1,4,6,21,24",
                "metadata_tsv_name": "sample_metadata_HOM_files_that_need_mixed_bcftools_concat_and_merge.tsv",
                "bucket_name": "mixed-concat-merge-project",
                "user_name": "test-user",
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
                "Cleaning up 10 temporary files",
            ],
            [],
        ),
        # case: expected to be sucessful, should lead to two different sample sets being concatenated, then merged with a third file
        (
            {
                "tsv_filter": "Area:Northern Spanish shelf,Iceland",
                "command": "view -s SAMPLES; view -r 1,4,6,8,13,18,21,24",
                "metadata_tsv_name": "sample_metadata_HOM_files_that_need_mixed_bcftools_concat_and_merge.tsv",
                "bucket_name": "mixed-concat-merge-project",
                "user_name": "test-user",
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
                "Cleaning up 17 temporary files",
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

        logger = logging.getLogger("divbase_worker.tasks")

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

        s3_client = boto3.client(
            "s3",
            endpoint_url=CONSTANTS["MINIO_URL"],
            aws_access_key_id=CONSTANTS["BAD_ACCESS_KEY"],
            aws_secret_access_key=CONSTANTS["BAD_SECRET_KEY"],
        )
        obj = s3_client.get_object(Bucket=params["bucket_name"], Key=DIMENSIONS_FILE_NAME)
        dimensions = yaml.safe_load(obj["Body"])
        print("VCF Dimensions file for bucket", params["bucket_name"])
        for filename, samples in dimensions.items():
            print(f"{filename}: {samples}")

    try:
        current_app.conf.update(
            task_always_eager=True,
            task_eager_propagates=True,
        )

        with (
            patch("boto3.client", side_effect=patched_boto3_client),
            patch("divbase_worker.tasks.download_sample_metadata", new=patched_download_sample_metadata),
            patch("divbase_lib.queries.BcftoolsQueryManager.CONTAINER_NAME", "divbase-tests-worker-quick-1"),
            patch("divbase_worker.tasks.download_vcf_files", side_effect=patched_download_vcf_files),
            patch("divbase_lib.queries.BcftoolsQueryManager.run_bcftools", new=patched_run_bcftools),
            patch("divbase_lib.queries.BcftoolsQueryManager.temp_file_management", new=patched_temp_file_management),
            patch(
                "divbase_lib.queries.BcftoolsQueryManager.merge_or_concat_bcftools_temp_files",
                new=patched_merge_or_concat_bcftools_temp_files,
            ),
            patch("divbase_worker.tasks.upload_results_file", new=patched_upload_results_file),
            patch("divbase_worker.tasks.delete_job_files_from_worker", new=patched_delete_job_files_from_worker),
        ):
            if not expect_success:
                with pytest.raises(ValueError) as e:
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


@patch(
    "divbase_worker.tasks.create_s3_file_manager",
    side_effect=lambda url=None: create_s3_file_manager(url="http://localhost:9002"),
)
def test_query_exits_when_vcf_file_version_is_outdated(
    mock_create_s3_manager,
    CONSTANTS,
    user_config_path,
    fixtures_dir,
    run_update_dimensions,
):
    """
    Test that updates the dimensions file, uploads a new version of a VCF file, then runs a query that should fail
    because the dimensions file expects an older version of the VCF file.
    """
    bucket_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    # this factory fixture needs to be run before the with patch below since otherwise the patfch will also affect the fixture
    run_update_dimensions(bucket_name=bucket_name)

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
        patch("divbase_worker.tasks.download_sample_metadata", new=patched_download_sample_metadata),
        patch("divbase_worker.tasks.download_vcf_files", new=patched_download_vcf_files),
    ):
        test_file = (fixtures_dir / "HOM_20ind_17SNPs.1.vcf.gz").resolve()

        command = f"files upload {test_file}  --project {bucket_name}"
        result = runner.invoke(app, command)

        assert result.exit_code == 0
        assert f"{str(test_file)}" in result.stdout

        params = {
            "tsv_filter": "Area:West of Ireland;Sex:F",
            "command": "view -s SAMPLES; view -r 1,4,6,21,24",
            "metadata_tsv_name": "sample_metadata_HOM_chr_split_version.tsv",
            "bucket_name": "split-scaffold-project",
            "user_name": "test-user",
        }
        with pytest.raises(ValueError) as excinfo:
            bcftools_pipe_task(**params)
        assert (
            "The VCF dimensions file is not up to date with the VCF files in the project. Please run 'divbase-cli dimensions update --project <project_name>' and then submit the query again."
            in str(excinfo.value)
        )
