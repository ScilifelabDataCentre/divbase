import datetime
import logging
import os
import subprocess
from contextlib import contextmanager
from pathlib import Path
from unittest.mock import patch

import pytest
import structlog
from celery import current_app
from kombu.connection import Connection
from structlog.testing import capture_logs
from typer.testing import CliRunner

from divbase_api.exceptions import VCFDimensionsEntryMissingError
from divbase_api.services.bcftools_helpers import get_container_id
from divbase_api.services.vcf_queries import BcftoolsQueryManager
from divbase_api.worker.tasks import bcftools_pipe_task
from divbase_cli.divbase_cli import app
from divbase_lib.divbase_constants import QUERY_RESULTS_FILE_PREFIX
from divbase_lib.exceptions import DimensionsNotUpToDateWithBucketError, TaskUserError
from tests.conftest import _text_in_logs

runner = CliRunner()


@pytest.fixture(autouse=True, scope="function")
def auto_clean_dimensions_entries_for_all_projects(clean_all_projects_dimensions):
    """Enable auto-cleanup of dimensions entries for all tests in this test file."""
    yield


@pytest.mark.integration
def test_bcftools_pipe_task_with_real_worker(
    wait_for_celery_task_completion,
    bcftools_pipe_kwargs_fixture,
    run_update_dimensions,
    project_map,
    CONSTANTS,
):
    """
    Integration test in which bcftools_pipe_task is run with a real Celery worker.
    Runs locally using the docker compose testing stack defined and performs a real sidecar and bcftools
    query by loading VCF files from the tests/fixtures dir. It was designed for having RabbitMQ as the broker, PostgreSQL as the backend,
    and a custom Celery worker image that has bcftools installed.
    (this test does not download fixture from bucket, since that is handled by the CLI layer)

    Does not assert file download and uploads, since that is handled by tests in tests/cli_commands/test_query_cli.py.
    (but the current version of the task does I/O with files from the bucket, so it is tested indirectly)

    """

    project_name = CONSTANTS["QUERY_PROJECT"]
    project_id = project_map[project_name]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    run_update_dimensions(bucket_name=bucket_name, project_id=project_id, project_name=project_name)
    bcftools_pipe_kwargs_fixture["project_id"] = project_id

    broker_url = current_app.conf.broker_url
    with Connection(broker_url) as conn:
        conn.ensure_connection(max_retries=1)

    async_result = bcftools_pipe_task.apply_async(kwargs=bcftools_pipe_kwargs_fixture)
    task_id = async_result.id
    task_result = wait_for_celery_task_completion(task_id=task_id, max_wait=30)

    assert task_result["status"] == "completed"
    assert task_result["output_file"].startswith(QUERY_RESULTS_FILE_PREFIX) and task_result["output_file"].endswith(
        ".vcf.gz"
    )


@pytest.mark.integration
@pytest.mark.parametrize(
    "tsv_filter,command,expected_error",
    [
        (
            "Area West of Ireland",
            "view -r 21:15000000-25000000",
            "Invalid filter format",
        ),
    ],
)
def test_bcftools_pipe_task_errors_task_failure(
    wait_for_celery_task_completion,
    run_update_dimensions,
    project_map,
    CONSTANTS,
    tsv_filter,
    command,
    expected_error,
):
    """Test task-time failures for bcftools pipe task execute path."""
    project_name = CONSTANTS["QUERY_PROJECT"]
    project_id = project_map[project_name]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    run_update_dimensions(bucket_name=bucket_name, project_id=project_id, project_name=project_name)

    broker_url = current_app.conf.broker_url
    with Connection(broker_url) as conn:
        conn.ensure_connection(max_retries=1)

    params = {
        "tsv_filter": tsv_filter,
        "command": command,
        "metadata_tsv_name": "sample_metadata.tsv",
        "bucket_name": bucket_name,
        "project_id": project_id,
        "project_name": project_name,
        "user_id": 1,
        "job_id": 19001,
    }
    async_result = bcftools_pipe_task.apply_async(kwargs=params)

    with pytest.raises(Exception) as exc_info:
        wait_for_celery_task_completion(task_id=async_result.id, max_wait=30)

    assert expected_error in str(exc_info.value)


@pytest.mark.parametrize(
    "test_scenario,vcf_filename,job_id,cleanup_file",
    [
        ("version_outdated", "HOM_20ind_17SNPs.1.vcf.gz", 1, False),
        ("unindexed", "HOM_20ind_17SNPs_first_10_samples.vcf.gz", 2, True),
    ],
)
def test_query_exits_when_dimensions_are_outdated(
    CONSTANTS,
    logged_in_edit_user_with_existing_config,
    fixtures_dir,
    run_update_dimensions,
    project_map,
    test_scenario,
    vcf_filename,
    job_id,
    cleanup_file,
):
    """
    Test that verifies DimensionsNotUpToDateWithBucketError is raised when the dimensions cache is not up-to-date. Test for these cases:
    1. version_outdated: uploads a new version of an existing VCF file after dimensions update
    2. unindexed: uploads a new VCF file that is not present in the dimensions cache
    """
    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    project_id = project_map[project_name]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    user_id = 1
    run_update_dimensions(bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id)

    def ensure_fixture_path(filename, fixture_dir="tests/fixtures"):
        if filename.startswith(fixture_dir):
            return filename
        return f"{fixture_dir}/{filename}"

    def patched_download_sample_metadata(metadata_tsv_name, bucket_name, s3_file_manager):
        return Path(ensure_fixture_path(metadata_tsv_name, fixture_dir="tests/fixtures"))

    def patched_download_vcf_files(files_to_download, bucket_name, s3_file_manager):
        pass

    with (
        patch("divbase_api.worker.tasks._download_sample_metadata", new=patched_download_sample_metadata),
        patch("divbase_api.worker.tasks._download_vcf_files", new=patched_download_vcf_files),
    ):
        test_file_path = (fixtures_dir / vcf_filename).resolve()
        command = f"files upload {test_file_path}  --project {project_name} --disable-safe-mode"
        result = runner.invoke(app, command)

        try:
            assert result.exit_code == 0
            assert vcf_filename in result.stdout

            params = {
                "tsv_filter": "Area:West of Ireland;Sex:F",
                "command": "view -r 1,4,6,21,24",
                "metadata_tsv_name": "sample_metadata_HOM_chr_split_version.tsv",
                "bucket_name": bucket_name,
                "project_id": project_id,
                "project_name": project_name,
                "user_id": 1,
                "job_id": job_id,
            }
            with pytest.raises(DimensionsNotUpToDateWithBucketError) as excinfo:
                bcftools_pipe_task(**params)
            assert (
                "The following VCF files or file versions in the project are not part of the project's VCF dimensions"
                in str(excinfo.value)
            )
        finally:
            # Ensure cleanup always runs for the unindexed scenario to avoid test pollution.
            if cleanup_file:
                command = f"files rm {vcf_filename}  --project {project_name}"
                runner.invoke(app, command)


@pytest.mark.integration
@pytest.mark.parametrize(
    "params,expect_success,ensure_dimensions_file,expected_logs,expected_error_msgs",
    [
        # Case 0: expected to be fail, vcf dimensions file is empty so the early check for the file in the task raises an error
        (
            {
                "tsv_filter": "Area:West of Ireland;Sex:F",
                "command": "view -r 1,4,6,21,24",
                "metadata_tsv_name": "sample_metadata_HOM_chr_split_version.tsv",
                "project_name": "split-scaffold-project",
                "user_id": 1,
            },
            False,
            False,
            ["Starting bcftools_pipe_task"],
            [
                "The VCF dimensions cache in project 'split-scaffold-project' is missing or empty. Please ensure that there are VCF files in the project and run:'divbase-cli dimensions update --project <project_name>'"
            ],
        ),
        # Case 1: expected to be sucessful, should lead to concat
        (
            {
                "tsv_filter": "Area:West of Ireland;Sex:F",
                "command": "view -r 1,4,6,21,24",
                "metadata_tsv_name": "sample_metadata_HOM_chr_split_version.tsv",
                "project_name": "split-scaffold-project",
                "user_id": 1,
            },
            True,
            True,
            [
                "Starting bcftools_pipe_task",
                "No unsupported sample sets found. Proceeding with bcftools pipeline.",
                "Sample names overlap between some temp files, will concat overlapping sets, then merge if needed and possible.",
                "Only one file remained after concatenation, renamed this file to",
                f"Sorting the results file to ensure proper order of variants. Final results are in '{QUERY_RESULTS_FILE_PREFIX}",
                "bcftools processing completed successfully",
            ],
            [],
        ),
        # Case 2: expected to fail since there are no scaffolds in the vcf files that match the -r query
        (
            {
                "tsv_filter": "Area:West of Ireland;Sex:F",
                "command": "view -r 31,34,36,321,324",
                "metadata_tsv_name": "sample_metadata_HOM_chr_split_version.tsv",
                "project_name": "split-scaffold-project",
                # project_id is added dynamically in the tests
                "user_id": 1,
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
                "command": "view -r 1,4,6,21,24",
                "metadata_tsv_name": "sample_metadata_HOM_chr_split_version.tsv",
                "project_name": "split-scaffold-project",
                # project_id is added dynamically in the tests
                "user_id": 1,
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
                f"Sorting the results file to ensure proper order of variants. Final results are in '{QUERY_RESULTS_FILE_PREFIX}",
                "bcftools processing completed successfully",
            ],
            [],
        ),
        # Case 4: expected to be sucessful, should lead to merge
        (
            {
                "tsv_filter": "Area:Northern Portugal",
                "command": "view -r 21:15000000-25000000",
                "metadata_tsv_name": "sample_metadata.tsv",
                "project_name": "query-project",
                # project_id is added dynamically in the tests
                "user_id": 1,
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
                f"Sorting the results file to ensure proper order of variants. Final results are in '{QUERY_RESULTS_FILE_PREFIX}",
                "bcftools processing completed successfully",
            ],
            [],
        ),
        # Case 5: expected to be sucessful, should lead one file being subset and renamed rather than bcftools merge/concat
        (
            {
                "tsv_filter": "Area:Northern Spanish shelf",
                "command": "view -r 1,4,6,21,24",
                "metadata_tsv_name": "sample_metadata_HOM_files_that_need_mixed_bcftools_concat_and_merge.tsv",
                "project_name": "mixed-concat-merge-project",
                # project_id is added dynamically in the tests
                "user_id": 1,
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
            ],
            [],
        ),
        # Case 6: expected to be sucessful, should lead to two different sample sets being concatenated, then merged with a third file
        (
            {
                "tsv_filter": "Area:Northern Spanish shelf,Iceland",
                "command": "view -r 1,4,6,8,13,18,21,24",
                "metadata_tsv_name": "sample_metadata_HOM_files_that_need_mixed_bcftools_concat_and_merge.tsv",
                "project_name": "mixed-concat-merge-project",
                # project_id is added dynamically in the tests
                "user_id": 1,
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
):
    """
    This is a special integration test that allows for running vcf queries directly in eager mode
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
    # Add job_id, incrementing for each parameterization
    if not hasattr(test_bcftools_pipe_cli_integration_with_eager_mode, "call_count"):
        test_bcftools_pipe_cli_integration_with_eager_mode.call_count = 1
    params["job_id"] = test_bcftools_pipe_cli_integration_with_eager_mode.call_count
    test_bcftools_pipe_cli_integration_with_eager_mode.call_count += 1

    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    caplog.set_level(logging.INFO)

    original_task_always_eager = current_app.conf.task_always_eager
    original_task_eager_propagates = current_app.conf.task_eager_propagates
    original_merge_or_concat_bcftools_temp_files = BcftoolsQueryManager.merge_or_concat_bcftools_temp_files

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
        return [
            Path(ensure_fixture_path(file_name, fixture_dir="/app/tests/fixtures")) for file_name in files_to_download
        ]

    def patched_run_bcftools(command: str, capture_output: bool = False, capture_stderr: bool = False):
        """
        Patches the working dir used when running bcftools commands inside the Docker container.
        """

        class DummyProc:
            """
            Dummy process object to mock subprocess.Popen for testing.
            Provides a .pid attribute and minimal methods so that code
            expecting a real process (for per-task resource monitoring/metrics)
            does not fail. Used to simulate a running process in tests
            where no actual subprocess is started.
            """

            def __init__(self, returncode: int = 0, stdout: str = "", stderr: str = ""):
                self.pid = 12345
                self.returncode = returncode
                self._stdout = stdout
                self._stderr = stderr

            # Simulate sucessful completion using 0
            def wait(self):
                return self.returncode

            def poll(self):
                return self.returncode

            def communicate(self):
                return self._stdout, self._stderr

        container_id = get_container_id("divbase-tests-worker-quick-1")
        logger = structlog.get_logger(__name__)
        logger.debug(f"Executing command in container with ID: {container_id}")
        docker_cmd = ["docker", "exec", "-w", "/app/tests/fixtures", container_id, "bcftools"] + command.split()
        run_result = subprocess.run(
            docker_cmd,
            check=not (capture_output or capture_stderr),
            capture_output=(capture_output or capture_stderr),
            text=True,
        )
        if capture_output or capture_stderr:
            return DummyProc(
                returncode=run_result.returncode,
                stdout=run_result.stdout or "",
                stderr=run_result.stderr or "",
            )
        return DummyProc()

    @contextmanager
    def patched_temp_file_management(self):
        """Context manager to handle temporary file cleanup, ensuring all temp files are in ./tests/fixtures."""
        self.temp_files = []
        try:
            yield self
        finally:
            if self.temp_files:
                logger = structlog.get_logger(__name__)
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
        original_group_vcfs_by_sample_set = self._group_vcfs_by_sample_set

        def patched_rename(src, dst):
            src = ensure_fixture_path(src)
            dst = ensure_fixture_path(dst)
            return original_rename(src, dst)

        def patched_get_all_sample_names_from_vcf_files(output_temp_files):
            """
            Extract sample names from temp VCF files by running bcftools in the worker container.
            This avoids relying on a host-installed bcftools binary when running eager-mode tests.
            """
            sample_names_per_VCF = {}
            container_id = get_container_id("divbase-tests-worker-quick-1")

            for vcf_file in output_temp_files:
                local_fixture_path = ensure_fixture_path(vcf_file)
                container_fixture_path = f"/app/tests/fixtures/{Path(local_fixture_path).name}"

                result = subprocess.run(
                    [
                        "docker",
                        "exec",
                        "-w",
                        "/app/tests/fixtures",
                        container_id,
                        "bcftools",
                        "query",
                        "-l",
                        container_fixture_path,
                    ],
                    capture_output=True,
                    text=True,
                    check=True,
                )
                sample_names = result.stdout.strip().split("\n")
                sample_names_per_VCF[local_fixture_path] = sample_names if sample_names != [""] else []

            return sample_names_per_VCF

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

    def patched_delete_job_files_from_worker(vcf_paths=None, metadata_path=None, output_file=None, log_file=None):
        """
        Only delete the output file, using the correct path. Don't delete the fixtures, since they should persist.
        """

        logger = structlog.get_logger(__name__)

        if output_file is not None:
            output_file = ensure_fixture_path(str(output_file))
            try:
                os.remove(output_file)
                logger.info(f"deleted {output_file}")
            except Exception as e:
                logger.warning(f"Could not delete output file from worker {output_file}: {e}")

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
            patch("divbase_api.worker.tasks._download_sample_metadata", new=patched_download_sample_metadata),
            patch("divbase_api.worker.tasks._download_vcf_files", side_effect=patched_download_vcf_files),
            patch("divbase_api.services.vcf_queries.run_bcftools", new=patched_run_bcftools),
            patch("divbase_api.services.bcftools_helpers.run_bcftools", new=patched_run_bcftools),
            patch(
                "divbase_api.services.vcf_queries.BcftoolsQueryManager.temp_file_management",
                new=patched_temp_file_management,
            ),
            patch(
                "divbase_api.services.vcf_queries.BcftoolsQueryManager.merge_or_concat_bcftools_temp_files",
                new=patched_merge_or_concat_bcftools_temp_files,
            ),
            patch("divbase_api.worker.tasks._upload_results_file", new=patched_upload_results_file),
            patch("divbase_api.worker.tasks._delete_job_files_from_worker", new=patched_delete_job_files_from_worker),
            patch(
                "divbase_api.services.vcf_queries.BcftoolsQueryManager._prepare_txt_with_divbase_header_for_vcf",
                new=patched_prepare_txt_with_divbase_header_for_vcf,
            ),
        ):
            with capture_logs() as cap_logs:
                if not expect_success:
                    with pytest.raises((VCFDimensionsEntryMissingError, TaskUserError)) as e:
                        bcftools_pipe_task(**params)
                    for msg in expected_error_msgs:
                        assert msg.replace("\n", "") in str(e.value).replace("\n", "")
                else:
                    result = bcftools_pipe_task(**params)
                    assert result is not None

                for log_msg in expected_logs:
                    assert _text_in_logs(text=log_msg, logs=cap_logs)

            print(f"Captured logs:\n{cap_logs}")
    finally:
        current_app.conf.task_always_eager = original_task_always_eager
        current_app.conf.task_eager_propagates = original_task_eager_propagates
