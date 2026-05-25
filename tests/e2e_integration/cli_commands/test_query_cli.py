"""
e2e tests for the "divbase-cli query" commands

All tests are run against a docker compose setup with the entire DivBase stack running locally.

A project (CONSTANTS["QUERY_PROJECT"]) is made available with input files for the tests.
"""

import csv
import gzip
import io
import re
import time
from pathlib import Path

import boto3
import pytest
from typer.testing import CliRunner

from divbase_cli.cli_exceptions import ProjectNotInConfigError
from divbase_cli.divbase_cli import app
from divbase_lib.divbase_constants import QUERY_RESULTS_FILE_PREFIX
from divbase_lib.s3_checksums import MD5CheckSumFormat, calculate_md5_checksum
from tests.conftest import REGRESSION_GUARD_PREFIX

runner = CliRunner()

_TASK_TERMINAL_STATES = {"SUCCESS", "FAILURE"}


@pytest.fixture(autouse=True, scope="function")
def auto_clean_dimensions_entries_for_all_projects(clean_all_projects_dimensions):
    """Enable auto-cleanup of dimensions entries for all tests in this test file."""
    yield


def _extract_results_file_key(text: str) -> str | None:
    """Extract the query results `.vcf.gz` key from user-visible CLI output."""
    match = re.search(rf"({re.escape(QUERY_RESULTS_FILE_PREFIX)}[^\s\"']*\.vcf\.gz)", text)
    return match.group(1) if match else None


def _extract_task_state_from_task_history_tsv(stdout: str) -> str | None:
    """
    Infer task state from `divbase-cli task-history id --tsv` output.

    Returns the state string from the TSV `State` column when present; otherwise None.
    """
    try:
        reader = csv.DictReader(io.StringIO(stdout), delimiter="\t")
        first_row = next(reader, None)
    except csv.Error:
        return None

    if first_row is None:
        return None

    state = first_row.get("State")
    if state is None:
        return None

    return state.strip().upper()


def wait_for_task_terminal_state_using_CLI(
    user_task_id: int | str, max_retries: int = 60, retry_delay: int = 1
) -> tuple[str, str]:
    """
    Repeatedly poll `divbase-cli task-history id --tsv` until the Celery task reaches one of the two terminal state
    used for DivBase celery tasks: SUCCESS or FAILURE.

    By using the CLI command to poll the task state instead of direct db lookups, this becomes a e2e test helper.
    """
    latest_stdout = ""
    transient_visibility_error = "Task ID not found or you don't have permission to view the history for this task ID."
    # For testing cases, we know that the task ID will exist, so if the signal hits this message, it just means to wait and try the CLI cmd again

    for attempt in range(max_retries):
        result = runner.invoke(app, f"task-history id {user_task_id} --tsv")
        latest_stdout = result.stdout

        if result.exit_code == 0:
            state = _extract_task_state_from_task_history_tsv(latest_stdout)
            if state in _TASK_TERMINAL_STATES:
                return state, latest_stdout
            time.sleep(retry_delay)
            continue

        # Right after submission, task-history can transiently return this before task metadata is queryable.
        combined_output = f"{latest_stdout}\n{result.exception!r}" if result.exception else latest_stdout
        if transient_visibility_error in combined_output:
            time.sleep(retry_delay)
            continue

        # Fail fast on unexpected CLI errors.
        raise AssertionError(
            "task-history polling failed unexpectedly.\n"
            f"task_id={user_task_id}, attempt={attempt + 1}/{max_retries}, exit_code={result.exit_code}\n"
            f"stdout:\n{latest_stdout}\n"
            f"exception: {result.exception!r}"
        )

    raise TimeoutError(
        f"Task {user_task_id} did not reach SUCCESS/FAILURE after {max_retries} retries. "
        f"Latest output:\n{latest_stdout}"
    )


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


def _checksum_vcf_skip_double_hash_headers(vcf_gz_file: Path, tmp_path: Path) -> str:
    """
    Calculate MD5 for a VCF file after stripping all `##` header lines. Intended for toy VCF results files since this implementation
    was not designed to be performant for scaling to large VCF files.

    The function requires a VCF file and not streamed bytes, since calculate_md5_checksum() expects a file path.

    VCF files contains ## and # headers. The challenge with the ## headers are that they can contain timestamped lines
    that change the checksum. This is part of the bcftools audit trail. It can be turned off, but that is not desired for DivBase.
    Timestamped lines will be part of the DivBase results files since the queries run 'bcftools view'. And since many tests
    will make assertions on the results files, we need to skip the timestamped lines if we want to make assertions on the checksums.
    Otherwise the checksums will not be deterministic.

    Example: the VCF file fixture in tests/fixtures/HOM_20ind_17SNPs_first_10_samples.vcf.gz contains the ## header
    ##bcftools_viewCommand=view -S subset_samples.txt -Oz -o HOM_20ind_17SNPs.vcf.gz ...Date=Tue May 27 18:16:01 2025

    There is also one # header line. That one contains the canonical columns and all the sample name columns. That one should be kept in
    for the checksum calculation.
    """
    text = gzip.decompress(vcf_gz_file.read_bytes()).decode("utf-8")
    kept_lines = [line for line in text.splitlines() if not line.startswith("##")]
    normalized_bytes = ("\n".join(kept_lines) + "\n").encode("utf-8")

    normalized_file = tmp_path / "vcf_skip_double_hash_headers.vcf"
    normalized_file.write_bytes(normalized_bytes)

    return calculate_md5_checksum(
        file_path=normalized_file,
        output_format=MD5CheckSumFormat.HEX,
    )


class TestQueryTSVSuccess:
    """Successful query tsv CLI test cases."""

    def test_sample_metadata_query(
        self,
        CONSTANTS,
        logged_in_edit_user_with_existing_config,
        run_update_dimensions,
        project_map,
    ):
        """Test running a sample metadata query using the CLI."""
        project_name = CONSTANTS["QUERY_PROJECT"]
        project_id = project_map[project_name]
        bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
        user_id = 1
        run_update_dimensions(
            bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
        )

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

    def test_sample_metadata_query_prints_explicit_message_when_no_samples_match(
        self,
        CONSTANTS,
        logged_in_edit_user_with_existing_config,
        run_update_dimensions,
        project_map,
    ):
        """Test that query tsv prints an explicit message when filters match zero samples."""
        project_name = CONSTANTS["QUERY_PROJECT"]
        project_id = project_map[project_name]
        bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
        user_id = 1
        run_update_dimensions(
            bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
        )

        query_string = "Area:ValueThatDoesNotExistInFixtures"
        command = f"query tsv '{query_string}' --project {project_name}"
        result = runner.invoke(app, command)

        assert result.exit_code == 0
        assert "Unique Sample IDs: []" in result.stdout
        assert "Unique filenames: []" in result.stdout
        assert "No samples match your query filters." in result.stdout


class TestQueryVCFSuccess:
    """Successful query vcf CLI test cases."""

    def test_bcftools_pipe_query(
        self,
        CONSTANTS,
        logged_in_edit_user_with_existing_config,
        run_update_dimensions,
        project_map,
    ):
        """Test running a bcftools pipe query using the CLI."""
        project_name = CONSTANTS["QUERY_PROJECT"]
        project_id = project_map[project_name]
        bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
        user_id = 1
        run_update_dimensions(
            bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
        )
        tsv_filter = "Area:West of Ireland,Northern Portugal;"
        arg_command = "view -r 21:15000000-25000000"

        command = f"query vcf --tsv-filter '{tsv_filter}' --command '{arg_command}' --project {project_name} "
        result = runner.invoke(app, command)
        assert result.exit_code == 0
        assert "Job submitted" in result.stdout

        user_task_id = result.stdout.strip().split()[-1]
        task_state, task_stdout = wait_for_task_terminal_state_using_CLI(user_task_id=user_task_id)
        assert task_state == "SUCCESS", f"Task failed. task-history output:\n{task_stdout}"

        command = f"files ls --project {project_name} --include-results-files"
        result = runner.invoke(app, command)

        assert result.exit_code == 0
        assert any(QUERY_RESULTS_FILE_PREFIX in line for line in result.stdout.splitlines()), (
            f"No {QUERY_RESULTS_FILE_PREFIX} VCF file found in output.\nfiles ls output:\n{result.stdout}"
        )

    @pytest.mark.parametrize(
        "files_to_upload, metadata_tsv_name, sample_selection_args, bcftools_view_command, expected_checksum",
        [
            (  # Case 1: Sample metadata HOM_20ind_17SNPs (the DivBase dev classic test case)
                [
                    "HOM_20ind_17SNPs_first_10_samples.vcf.gz",
                    "HOM_20ind_17SNPs_last_10_samples.vcf.gz",
                    "sample_metadata_HOM_20ind_17SNPs.tsv",
                ],
                "sample_metadata_HOM_20ind_17SNPs.tsv",
                "--tsv-filter 'Area:West of Ireland,Northern Portugal;Sex:F'",
                "view -s",
                "3f9c371bcffb8126663cf08a802ae58c",
            ),
            (  # Case 2: Sample metadata HOM_chr_split_version.
                # Interestingly, while the data content for this results filecase is the same of case 1, it uses bcftools concat since split VCF files and thus the sample column order is different => different checksum.
                [
                    "HOM_20ind_17SNPs.1.vcf.gz",
                    "HOM_20ind_17SNPs.4.vcf.gz",
                    "HOM_20ind_17SNPs.5.vcf.gz",
                    "HOM_20ind_17SNPs.6.vcf.gz",
                    "HOM_20ind_17SNPs.7.vcf.gz",
                    "HOM_20ind_17SNPs.8.vcf.gz",
                    "HOM_20ind_17SNPs.13.vcf.gz",
                    "HOM_20ind_17SNPs.18.vcf.gz",
                    "HOM_20ind_17SNPs.20.vcf.gz",
                    "HOM_20ind_17SNPs.21.vcf.gz",
                    "HOM_20ind_17SNPs.22.vcf.gz",
                    "HOM_20ind_17SNPs.24.vcf.gz",
                    "sample_metadata_HOM_20ind_17SNPs.tsv",
                ],
                "sample_metadata_HOM_20ind_17SNPs.tsv",
                "--tsv-filter 'Area:West of Ireland,Northern Portugal;Sex:F'",
                "view -s",
                "f44ecb1c1ebc5e03a0e9055a699b8cee",
            ),
            (  # Case 3: Sample metadata HOM_chr_split_version with view -r subset. Based on intergration test cases in tests/e2e_integration/queries_integration/test_bcftools_tasks.py.
                [
                    "HOM_20ind_17SNPs.1.vcf.gz",
                    "HOM_20ind_17SNPs.4.vcf.gz",
                    "HOM_20ind_17SNPs.5.vcf.gz",
                    "HOM_20ind_17SNPs.6.vcf.gz",
                    "HOM_20ind_17SNPs.7.vcf.gz",
                    "HOM_20ind_17SNPs.8.vcf.gz",
                    "HOM_20ind_17SNPs.13.vcf.gz",
                    "HOM_20ind_17SNPs.18.vcf.gz",
                    "HOM_20ind_17SNPs.20.vcf.gz",
                    "HOM_20ind_17SNPs.21.vcf.gz",
                    "HOM_20ind_17SNPs.22.vcf.gz",
                    "HOM_20ind_17SNPs.24.vcf.gz",
                    "sample_metadata_HOM_20ind_17SNPs.tsv",
                ],
                "sample_metadata_HOM_20ind_17SNPs.tsv",
                "--tsv-filter 'Area:West of Ireland;Sex:F'",
                "view -r 1,4,6,21,24",
                "9ef2d87a9adf577c2f45b47d1e70138b",
            ),
        ],
        ids=[
            "case-1-sample-metadata-HOM-20ind",
            "case-2-sample-metadata-HOM-chr-split",
            "case-3-sample-metadata-HOM-chr-split-view-r-subset",
        ],
    )
    def test_vcf_query_result_file_by_headerless_checksum(
        self,
        CONSTANTS,
        logged_in_edit_user_with_existing_config,
        run_update_dimensions,
        project_map,
        fixtures_dir,
        cleaned_project_bucket,
        files_to_upload,
        metadata_tsv_name,
        sample_selection_args,
        bcftools_view_command,
        expected_checksum,
        tmp_path,
    ):
        """
        Test that vcf queries for a given project, sample selection command, and bcftools view command produce the expected ##-headerless VCF checksum.

        The expected checksum should be calculated with _checksum_vcf_skip_double_hash_headers() to skip the ## headers (since they can contain timestamps that change the checksum).
        For more details on how to calculate the expected checksum, see docs/development/writing_e2e_tests_for_vcf_results_checksums.md.

        The parametrised test cases for this test might overlap with other e2e and integration test cases; this is intentional as this test aims to assert on stable results file output over time.
        Other tests might assert on more granular details of the results file generation.
        """
        # The cleaned_project_bucket fixture is function-scoped, so its setup/teardown runs once
        # per parametrized test invocation, keeping this test isolated.
        project_name, bucket_name = cleaned_project_bucket
        assert project_name == CONSTANTS["CLEANED_PROJECT"]

        project_id = project_map[project_name]
        user_id = 1

        for fixture_name in files_to_upload:
            fixture_path = (fixtures_dir / fixture_name).resolve()
            upload_result = runner.invoke(app, f"files upload {fixture_path} --project {project_name}")
            assert upload_result.exit_code == 0, f"Upload failed for {fixture_name}: {upload_result.stdout}"

        run_update_dimensions(
            bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
        )

        metadata_arg = f" --metadata-tsv-name {metadata_tsv_name}" if metadata_tsv_name else ""
        query_result = runner.invoke(
            app,
            f"query vcf {sample_selection_args} --command '{bcftools_view_command}' --project {project_name}{metadata_arg}",
        )
        assert query_result.exit_code == 0, f"Query submission failed: {query_result.stdout}"
        assert "Job submitted" in query_result.stdout

        user_task_id = query_result.stdout.strip().split()[-1]
        task_state, task_stdout = wait_for_task_terminal_state_using_CLI(user_task_id=user_task_id)
        assert task_state == "SUCCESS", f"Expected SUCCESS. task-history output:\n{task_stdout}"

        output_file = f"{QUERY_RESULTS_FILE_PREFIX}{user_task_id}.vcf.gz"

        # The _checksum_vcf_skip_double_hash_headers eventually requires a VCF file and not streamed bytes, so we can either download it or stream it to a file.
        # The benefit of using divbase-cli files download is that that command runs a download integrity checksum test, meaning that we know that the file that
        # will be stripped for ## headers by _checksum_vcf_skip_double_hash_headers is intact at the helper function call time.
        download_result = runner.invoke(
            app,
            f"files download {output_file} --download-dir {tmp_path} --project {project_name}",
        )
        assert download_result.exit_code == 0, f"Downloading query result failed: {download_result.stdout}"

        downloaded_results_file = tmp_path / output_file
        assert downloaded_results_file.exists(), f"Expected downloaded file not found: {downloaded_results_file}"

        checksum = _checksum_vcf_skip_double_hash_headers(vcf_gz_file=downloaded_results_file, tmp_path=tmp_path)
        assert checksum == expected_checksum, f"Expected checksum {expected_checksum}, got {checksum}"

    @pytest.mark.parametrize(
        "arg_command,expected_records,expected_view_command_filter_fragment",
        [
            (
                r"view -i 'FILTER~\"q10\"'",
                [("17330", ".", "q10"), ("1234567", "microsat1", "q10;s50")],
                'FILTER~"q10"',
            ),
            (
                r"view -i 'FILTER=\"q10\"'",
                [("17330", ".", "q10")],
                'FILTER="q10"',
            ),
            (
                r"view -i 'FILTER=\"q10;s50\"'",
                [("1234567", "microsat1", "q10;s50")],
                'FILTER="q10;s50"',
            ),
        ],
        ids=[
            "filter-subset-match-q10",
            "filter-exact-match-q10",
            "filter-exact-match-q10-s50",
        ],
    )
    def test_bcftools_pipe_query_supports_semicolon_in_filter_expression_e2e(
        self,
        CONSTANTS,
        logged_in_edit_user_with_existing_config,
        run_update_dimensions,
        project_map,
        fixtures_dir,
        cleaned_project_bucket,
        arg_command,
        expected_records,
        expected_view_command_filter_fragment,
    ):
        """
        Test that VCF queries can support bcftools view -i FILTER expressions (including semicolon values)

        The fixture in "vcf_specification_v45_example11.vcf.gz" was taken from https://samtools.github.io/hts-specs/VCFv4.5.pdf
        page, 4 1.1 An example. One line was modified to include the case of "q10;s50".
        """
        project_name, bucket_name = cleaned_project_bucket
        project_id = project_map[project_name]
        user_id = 1
        fixture_name = "vcf_specification_v45_example11.vcf.gz"
        fixture_path = (fixtures_dir / fixture_name).resolve()
        upload_result = runner.invoke(app, f"files upload {fixture_path} --project {project_name}")
        assert upload_result.exit_code == 0, f"Upload failed: {upload_result.stdout}"

        run_update_dimensions(
            bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
        )

        query_result = runner.invoke(
            app,
            f'query vcf --samples "NA00001,NA00002,NA00003" --command "{arg_command}" --project {project_name}',
        )
        assert query_result.exit_code == 0, f"Query submission failed: {query_result.stdout}"
        assert "Job submitted" in query_result.stdout

        user_task_id = query_result.stdout.strip().split()[-1]
        task_state, task_stdout = wait_for_task_terminal_state_using_CLI(user_task_id=user_task_id)
        assert task_state == "SUCCESS", f"Task failed. task-history output:\n{task_stdout}"

        output_file = _extract_results_file_key(task_stdout)

        stream_result = runner.invoke(
            app,
            f"files stream {output_file} --project {project_name}",
        )
        assert stream_result.exit_code == 0, f"Streaming query result failed: {stream_result.stdout}"

        streamed_vcf_content = gzip.decompress(stream_result.stdout_bytes).decode("utf-8")
        assert "##bcftools_viewCommand=" in streamed_vcf_content
        assert expected_view_command_filter_fragment in streamed_vcf_content

        records = [line for line in streamed_vcf_content.splitlines() if line and not line.startswith("#")]
        parsed_records = [(cols[1], cols[2], cols[6]) for cols in (record.split("\t") for record in records)]

        assert parsed_records == expected_records, (
            f"Unexpected records for command '{arg_command}'. Expected {expected_records}, got {parsed_records}"
        )

    def test_bcftools_pipe_query_direct_samples_mode(
        self,
        CONSTANTS,
        logged_in_edit_user_with_existing_config,
        run_update_dimensions,
        project_map,
    ):
        """Test running a bcftools pipe query using direct sample IDs from CLI."""
        project_name = CONSTANTS["QUERY_PROJECT"]
        project_id = project_map[project_name]
        bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
        user_id = 1
        run_update_dimensions(
            bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
        )

        arg_command = "view -r 21:15000000-25000000"
        command = f"query vcf --samples '5a_HOM-I13,1b_HOM-G58' --command '{arg_command}' --project {project_name} "
        result = runner.invoke(app, command)

        assert result.exit_code == 0
        assert "Job submitted" in result.stdout
        user_task_id = result.stdout.strip().split()[-1]
        task_state, task_stdout = wait_for_task_terminal_state_using_CLI(user_task_id=user_task_id)
        assert task_state == "SUCCESS", f"Expected SUCCESS. task-history output:\n{task_stdout}"

    def test_bcftools_pipe_query_all_samples_mode(
        self,
        CONSTANTS,
        logged_in_edit_user_with_existing_config,
        run_update_dimensions,
        project_map,
    ):
        """Test running a bcftools pipe query with explicit --all-samples mode."""
        project_name = CONSTANTS["QUERY_PROJECT"]
        project_id = project_map[project_name]
        bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
        user_id = 1
        run_update_dimensions(
            bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
        )

        arg_command = "view -r 21:15000000-25000000"
        command = f"query vcf --all-samples --command '{arg_command}' --project {project_name} "
        result = runner.invoke(app, command)

        assert result.exit_code == 0
        assert "Job submitted" in result.stdout
        user_task_id = result.stdout.strip().split()[-1]
        task_state, task_stdout = wait_for_task_terminal_state_using_CLI(user_task_id=user_task_id)
        assert task_state == "SUCCESS", f"Expected SUCCESS. task-history output:\n{task_stdout}"

    def test_bcftools_pipe_query_succeeds_twice_without_dimensions_update_between_runs(
        self,
        CONSTANTS,
        logged_in_edit_user_with_existing_config,
        run_update_dimensions,
        project_map,
    ):
        """
        Tests that running a vcf query twice works without running dimensions update between the two runs.
        Results files are prefixed with QUERY_RESULTS_FILE_PREFIX and _check_that_dimensions_is_up_to_date_with_VCF_files_in_bucket
        skips files that begin with that prefix.
        """
        project_name = CONSTANTS["QUERY_PROJECT"]
        project_id = project_map[project_name]
        bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
        user_id = 1
        run_update_dimensions(
            bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
        )

        tsv_filter = "Area:West of Ireland,Northern Portugal;"
        arg_command = "view -r 21:15000000-25000000"
        command = f"query vcf --tsv-filter '{tsv_filter}' --command '{arg_command}' --project {project_name} "

        first_result = runner.invoke(app, command)
        assert first_result.exit_code == 0, f"First run failed: {first_result.stdout}"
        first_task_id = first_result.stdout.strip().split()[-1]
        first_task_state, first_task_stdout = wait_for_task_terminal_state_using_CLI(user_task_id=first_task_id)
        assert first_task_state == "SUCCESS", f"First run task failed. task-history output:\n{first_task_stdout}"

        second_result = runner.invoke(app, command)
        assert second_result.exit_code == 0, f"Second run failed: {second_result.stdout}"
        second_task_id = second_result.stdout.strip().split()[-1]
        second_task_state, second_task_stdout = wait_for_task_terminal_state_using_CLI(user_task_id=second_task_id)
        assert second_task_state == "SUCCESS", (
            f"Second run task failed. task-history output:\n{second_task_stdout}\n"
            "This suggests that result files from the first run are incorrectly triggering DimensionsNotUpToDateWithBucketError."
        )


class TestQueryVCFSubmissionValidation:
    """Validation and submission-failure scenarios for query vcf CLI."""

    def test_bcftools_pipe_query_without_selection_mode_fails_on_submission(
        self,
        CONSTANTS,
        logged_in_edit_user_with_existing_config,
        run_update_dimensions,
        project_map,
    ):
        """Test that query vcf without any sample-selection mode fails before job submission."""
        project_name = CONSTANTS["QUERY_PROJECT"]
        project_id = project_map[project_name]
        bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
        user_id = 1
        run_update_dimensions(
            bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
        )

        command = f"query vcf --command 'view -r 21:15000000-25000000' --project {project_name} "
        result = runner.invoke(app, command)

        assert result.exit_code == 2  # Typer exits with code 2 for argument validation errors
        assert isinstance(result.exception, SystemExit)
        output = result.stdout + (str(result.exception) if result.exception else "")
        assert "Job submitted successfully with task id:" not in output

    def test_bcftools_pipe_fails_on_project_not_in_config(self, CONSTANTS, logged_in_edit_user_with_existing_config):
        project_name = "non_existent_project"
        tsv_filter = "Area:West of Ireland,Northern Portugal;"
        arg_command = "view -r 21:15000000-25000000"

        command = f"query vcf --tsv-filter '{tsv_filter}' --command '{arg_command}' --project {project_name} "
        result = runner.invoke(app, command)
        assert isinstance(result.exception, ProjectNotInConfigError)

    @pytest.mark.parametrize(
        "project_name,tsv_filter,command,expected_error",
        [
            # bad command
            ("DEFAULT", "DEFAULT", "invalid-command", "Unsupported bcftools command"),
            # empty command string
            ("DEFAULT", "DEFAULT", "", "non-empty bcftools view string"),
            # sample-file options in --command are rejected early by task guard
            ("DEFAULT", "DEFAULT", "view -S samples.txt", "-S/--samples-file"),
            # output options in --command are rejected early by API guard
            ("DEFAULT", "DEFAULT", "view -Oz", "-O/--output-type"),
            # input VCF filenames in --command are rejected early by API guard
            (
                "DEFAULT",
                "DEFAULT",
                "view file.vcf.gz",
                "Do not provide VCF/BCF input filenames in '--command'",
            ),
            (
                "DEFAULT",
                "DEFAULT",
                "view -r file.vcf.gz",
                "Do not provide VCF/BCF input filenames in '--command'",
            ),
        ],
    )
    def test_bcftools_pipe_query_errors(
        self,
        run_update_dimensions,
        project_map,
        project_name,
        tsv_filter,
        command,
        expected_error,
        CONSTANTS,
        logged_in_edit_user_with_existing_config,
    ):
        """Test submission-time command validation errors for bcftools queries."""
        if "DEFAULT" in project_name:
            project_name = CONSTANTS["QUERY_PROJECT"]
        if "DEFAULT" in tsv_filter:
            tsv_filter = "Area:West of Ireland,Northern Portugal;"
        if "DEFAULT" in command:
            command = "view -r 21:15000000-25000000"

        project_id = project_map[project_name]
        bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
        user_id = 1
        run_update_dimensions(
            bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
        )

        command_str = f"query vcf --tsv-filter '{tsv_filter}' --command '{command}' --project {project_name} "

        response = runner.invoke(app, command_str)
        assert response.exit_code == 1

        output = response.stdout
        if response.exception is not None:
            output = f"{output}\n{response.exception}"
        normalized_output = " ".join(output.split())
        assert expected_error in normalized_output, (
            f"Expected '{expected_error}' in submission failure output, but got: {normalized_output}"
        )
        assert "Job submitted successfully with task id:" not in output

    def test_query_vcf_submission_prints_task_history_hint(
        self, CONSTANTS, logged_in_edit_user_with_existing_config, run_update_dimensions, project_map
    ):
        """Query submission should return a task id and a task-history follow-up command hint."""
        project_name = CONSTANTS["QUERY_PROJECT"]
        project_id = project_map[project_name]
        bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
        user_id = 1
        run_update_dimensions(
            bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
        )

        tsv_filter = "Area:West of Ireland,Northern Portugal;"
        arg_command = "view -r 21:15000000-25000000"

        command = f"query vcf --tsv-filter '{tsv_filter}' --command '{arg_command}' --project {project_name} "
        submit_result = runner.invoke(app, command)
        assert submit_result.exit_code == 0
        assert "Job submitted successfully with task id:" in submit_result.stdout
        assert "task-history id" in " ".join(submit_result.stdout.split())
        task_id = submit_result.stdout.strip().split()[-1]
        assert task_id.isdigit(), f"Expected numeric task id in submission output, got: {submit_result.stdout}"


class TestQueryTSVPreconditionErrors:
    """Precondition and setup-related errors for query tsv."""

    def test_error_in_terminal_for_sample_metadata_query_on_tsv_not_in_bucket(
        self,
        CONSTANTS,
        run_update_dimensions,
        project_map,
        logged_in_edit_user_with_existing_config,
    ):
        """
        Test that the sample metadata query raises the correct error when the specified TSV file is not found in the project bucket.

        Indirectly covers the case of the user misspells the TSV filename.
        """
        project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
        bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
        project_id = project_map[project_name]
        user_id = 1

        run_update_dimensions(
            bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
        )

        tsv_filename = "filename_that_does_not_exist_in_bucket.tsv"
        command = f'query tsv "Area:North" --metadata-tsv-name {tsv_filename} --project {project_name}'
        cli_result = runner.invoke(app, command)

        assert f"The sample metadata TSV file '{tsv_filename}' was not found in your project '{project_name}'" in str(
            cli_result.exception
        ), "Expected error message about missing TSV file in project bucket"

    def test_error_in_terminal_for_sample_metadata_query_when_no_dimensions_file(
        self,
        CONSTANTS,
        logged_in_edit_user_with_existing_config,
    ):
        """
        Test that the sample metadata query raises the correct error when there is no dimensions file in the project bucket.
        """
        project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
        # To test for this, do not run the update_dimensions fixture

        tsv_filename = "sample_metadata_HOM_chr_split_version.tsv"
        command = f'query tsv "Area:North" --metadata-tsv-name {tsv_filename} --project {project_name} '
        cli_result = runner.invoke(app, command)

        assert f"The VCF dimensions index in project '{project_name}' is missing or empty" in str(
            cli_result.exception
        ), "Expected error message about missing VCF dimensions file in project bucket"

    def test_error_in_terminal_for_sample_metadata_query_tsv_missing_should_be_raised_before_dimensions_check(
        self,
        CONSTANTS,
        logged_in_edit_user_with_existing_config,
    ):
        """
        Test that the missing TSV in bucket error is raised before the dimensions file check for a case when both are incorrect.
        """
        project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
        # To test for this, do not run the update_dimensions fixture

        filename = "sample_metadata.tsv"  # not in the split-scaffold-project bucket
        command = f'query tsv "Area:North" --metadata-tsv-name {filename} --project {project_name} '
        cli_result = runner.invoke(app, command)

        assert f"The sample metadata TSV file '{filename}' was not found in your project '{project_name}'" in str(
            cli_result.exception
        ), "Expected error message about missing TSV file in project bucket"


class TestQueryTSVContentValidation:
    """TSV content/syntax validation behavior for query tsv."""

    def test_error_in_terminal_for_invalid_filter_syntax(
        self,
        CONSTANTS,
        run_update_dimensions,
        project_map,
        logged_in_edit_user_with_existing_config,
    ):
        """
        Test that SidecarInvalidFilterError is raised and propagated to terminal when filter syntax is invalid.
        """
        project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
        bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
        project_id = project_map[project_name]
        user_id = 1

        run_update_dimensions(
            bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
        )

        tsv_filename = "sample_metadata_HOM_chr_split_version.tsv"

        invalid_filter = "Area North"  # Use invalid filter syntax (missing colon)
        command = f'query tsv "{invalid_filter}" --metadata-tsv-name {tsv_filename} --project {project_name}'
        cli_result = runner.invoke(app, command)

        assert f"Invalid filter format: '{invalid_filter}'. Expected format 'key:value1,value2' or" in str(
            cli_result.exception
        ), "Expected error message about invalid filter syntax"

    def test_error_in_terminal_when_querying_nonexistent_column_in_tsv(
        self,
        CONSTANTS,
        run_update_dimensions,
        project_map,
        logged_in_edit_user_with_existing_config,
    ):
        """
        Test that SidecarColumnNotFoundError is raised and propagated to terminal when querying non-existent column.
        """
        project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
        bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
        project_id = project_map[project_name]
        user_id = 1

        run_update_dimensions(
            bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
        )

        tsv_filename = "sample_metadata_HOM_chr_split_version.tsv"

        invalid_filter = "NonExistentColumn:value"  # Query a column that doesn't exist in the TSV
        command = f'query tsv "{invalid_filter}" --metadata-tsv-name {tsv_filename} --project {project_name}'
        cli_result = runner.invoke(app, command)

        output = cli_result.stdout + (str(cli_result.exception) if cli_result.exception else "")
        # Normalize whitespace to handle line wrapping
        normalized_output = " ".join(output.split())
        assert cli_result.exit_code == 1, "Expected exit code 1 for invalid filter condition"
        assert (
            "Invalid filter conditions: no valid filter conditions could be parsed from 'NonExistentColumn:value'"
            in normalized_output
        )
        assert "NonExistentColumn:value" in normalized_output, "Expected filter string in error message"

    def test_error_in_terminal_when_duplicate_sample_IDs_in_tsv(
        self,
        CONSTANTS,
        run_update_dimensions,
        project_map,
        logged_in_edit_user_with_existing_config,
        tmp_path,
    ):
        """
        Test that SidecarSampleIDError is raised and propagated to terminal when TSV has duplicate Sample_IDs.
        """
        project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
        bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
        project_id = project_map[project_name]
        user_id = 1

        run_update_dimensions(
            bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
        )

        tsv_file = tmp_path / "test_duplicate_sample_ids.tsv"
        tsv_file.write_text("#Sample_ID\tArea\nS1\tNorth\nS1\tSouth\n8_HOM-E59\tEast\n")
        command = f"files upload {tsv_file} --project {project_name}"
        result = runner.invoke(app, command)
        assert result.exit_code == 0

        command = f'query tsv "Area:North" --metadata-tsv-name {tsv_file.name} --project {project_name}'
        cli_result = runner.invoke(app, command)

        error_text = str(cli_result.exception)
        assert "Duplicate Sample_IDs found" in error_text and "Each Sample_ID must be unique." in error_text, (
            f"Expected error message about duplicate Sample_IDs, got: {error_text}"
        )

    def test_error_in_terminal_for_comma_in_metadata(
        self,
        CONSTANTS,
        run_update_dimensions,
        project_map,
        logged_in_edit_user_with_existing_config,
        tmp_path,
    ):
        """
        Test that TSV files with commas generate warnings (not errors) and can still be queried.
        Columns with commas are treated as string columns.
        """
        project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
        bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
        project_id = project_map[project_name]
        user_id = 1

        run_update_dimensions(
            bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
        )

        tsv_file = tmp_path / "test_comma_in_data.tsv"
        tsv_file.write_text("#Sample_ID\tPopulation\n8_HOM-E57\t1,2\n8_HOM-E59\t3\n")
        command = f"files upload {tsv_file} --project {project_name}"
        result = runner.invoke(app, command)
        assert result.exit_code == 0

        command = f'query tsv "Population:2" --metadata-tsv-name {tsv_file.name} --project {project_name}'
        cli_result = runner.invoke(app, command)

        assert cli_result.exit_code == 0, f"Query should succeed with comma warning. Output: {cli_result.output}"
        assert "comma" in cli_result.output.lower(), "Expected warning message about comma in metadata value"


def test_regression_vcf_query_fails_with_descriptive_error_for_malformed_sample_columns(
    CONSTANTS,
    logged_in_edit_user_with_existing_config,
    run_update_dimensions,
    project_map,
    fixtures_dir,
    cleaned_project_bucket,
):
    """
    Regression test (negative outcome): test that a VCF query run against a malformed VCF file where one variant row has
    fewer sample columns than the header, fails at task time with a descriptive TaskUserError and not just with
    BcftoolsCommandError("Process exited with code 1").

    Why: bcftools index is used a guard against several VCF formatting issues in DivBase, but this
    particular malformation is not caught by bcftools index. Until another VCF validation system is
    implemented in DivBase, this test checks that the error is user-friendly and propagated.

    Reference: docs/development/bcftools_task_constraints.md ("1.4. What is required in the header?").
    The bcftools stderr (exit code 8):
    [E::vcf_parse_format_check7] Number of columns at 20:14370 does not match the number of samples (2 vs 3)
    Error: VCF parse error
    """

    project_name, bucket_name = cleaned_project_bucket
    assert project_name == CONSTANTS["CLEANED_PROJECT"]
    project_id = project_map[project_name]
    user_id = 1

    fixture_name = "vcf_specification_v45_example11_incorrect_sample_column_count_mismatch.vcf.gz"
    fixture_path = (fixtures_dir / fixture_name).resolve()
    assert fixture_path.exists(), f"Missing fixture file: {fixture_path}"

    upload_result = runner.invoke(app, f"files upload {fixture_path} --project {project_name}")
    assert upload_result.exit_code == 0, f"Upload failed: {upload_result.stdout}"

    run_update_dimensions(
        bucket_name=bucket_name,
        project_id=project_id,
        project_name=project_name,
        user_id=user_id,
    )

    query_result = runner.invoke(
        app,
        f'query vcf --samples "NA00001,NA00002,NA00003" --command "view -s" --project {project_name}',
    )
    assert query_result.exit_code == 0, (
        f"{REGRESSION_GUARD_PREFIX} query submission should succeed even for a malformed VCF "
        f"(column mismatches are not detectable before task execution): {query_result.stdout}"
    )
    assert "Job submitted" in query_result.stdout

    user_task_id = query_result.stdout.strip().split()[-1]
    task_state, task_stdout = wait_for_task_terminal_state_using_CLI(user_task_id=user_task_id)

    assert task_state == "FAILURE", (
        f"{REGRESSION_GUARD_PREFIX} expected task FAILURE for malformed VCF column mismatch, got {task_state}. "
        f"task-history output:\n{task_stdout}"
    )
    assert "malformed VCF record lines" in task_stdout, (
        f"{REGRESSION_GUARD_PREFIX} expected descriptive 'malformed VCF record lines' message in task-history "
        f"instead of the opaque 'Process exited with code' error. task-history output:\n{task_stdout}"
    )


def test_regression_query_fails_for_vcf_with_partly_overlapping_sample_sets(
    CONSTANTS,
    logged_in_edit_user_with_existing_config,
    run_update_dimensions,
    project_map,
    fixtures_dir,
    cleaned_project_bucket,
):
    """
    Regression test (negative outcome): test that a VCF query against two files with partially overlapping
    sample sets fail with a descriptive TaskUserError before any VCF files are downloaded.

    Why: bcftools merge does not accept partially overlapping samples between VCF files. DivBase deliberately
    omits --force-samples from its bcftools merge calls, since that could result in undesired renaming of
    conflicting samples (e.g. S3 → 2:S3).

    Reference: docs/development/bcftools_task_constraints.md ("3.2. Sample names must be unique
    across all input files" and "2. How DivBase chooses between bcftools merge and concat").
    """

    project_name, bucket_name = cleaned_project_bucket
    assert project_name == CONSTANTS["CLEANED_PROJECT"]
    project_id = project_map[project_name]
    user_id = 1

    # Both files contain '5a_HOM-I7' but have different remaining samples — a partial overlap.
    fixture_names = [
        "HOM_20ind_17SNPs_first_10_samples.vcf.gz",
        "HOM_20ind_17SNPs_last_10_samples_with_edit_to_break_merge.vcf.gz",
    ]
    for fixture_name in fixture_names:
        fixture_path = (fixtures_dir / fixture_name).resolve()
        upload_result = runner.invoke(app, f"files upload {fixture_path} --project {project_name}")
        assert upload_result.exit_code == 0, f"Upload failed for {fixture_name}: {upload_result.stdout}"

    run_update_dimensions(bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id)

    # '5a_HOM-I7' appears in both fixtures, so call a sample selection on that to trigger the overlap validation.
    query_result = runner.invoke(
        app,
        f'query vcf --samples "5a_HOM-I7" --command "view -s" --project {project_name}',
    )
    assert query_result.exit_code == 0, (
        f"{REGRESSION_GUARD_PREFIX} query submission should succeed even for incompatible sample sets "
        f"(incompatibility is only detectable once the task runs): {query_result.stdout}"
    )
    assert "Job submitted" in query_result.stdout

    user_task_id = query_result.stdout.strip().split()[-1]
    task_state, task_stdout = wait_for_task_terminal_state_using_CLI(user_task_id=user_task_id)

    assert task_state == "FAILURE", (
        f"{REGRESSION_GUARD_PREFIX} expected task FAILURE for partly overlapping sample sets, got {task_state}. "
        f"task-history output:\n{task_stdout}"
    )
    assert "Sample sets that are partly overlapping" in task_stdout, (
        f"{REGRESSION_GUARD_PREFIX} expected descriptive 'Sample sets that are partly overlapping' message "
        f"in task-history. task-history output:\n{task_stdout}"
    )
