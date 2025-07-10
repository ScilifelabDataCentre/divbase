"""
e2e tests for the "divbase-cli query" commands

All tests are run against a docker compose setup with the entire DivBase stack running locally.

A "query" bucket is made available with input files for the tests.
"""

import shlex
import time
from unittest.mock import patch

import boto3
import pytest
from typer.testing import CliRunner

from divbase_tools.divbase_cli import app

runner = CliRunner()


@pytest.fixture(autouse=True)
def use_test_api_url():
    with patch("divbase_tools.cli_commands.query_cli.DIVBASE_API_URL", "http://localhost:8001"):
        yield


@pytest.fixture(autouse=True)
def reset_query_bucket(CONSTANTS):
    """Delete the merged.vcf.gz file from the query bucket before each test."""
    s3_resource = boto3.resource(
        "s3",
        endpoint_url=CONSTANTS["MINIO_URL"],
        aws_access_key_id=CONSTANTS["BAD_ACCESS_KEY"],
        aws_secret_access_key=CONSTANTS["BAD_SECRET_KEY"],
    )

    # pylance does not understand boto3 resource returns types, hence ignore below
    bucket = s3_resource.Bucket(CONSTANTS["QUERY_BUCKET"])  # type: ignore
    bucket.object_versions.filter(Prefix="merged.vcf.gz").delete()
    yield


def test_sample_metadata_query(CONSTANTS):
    """Test running a sample metadata query using the CLI."""
    bucket = CONSTANTS["QUERY_BUCKET"]
    query_string = "Area:West of Ireland,Northern Portugal;Sex:F"
    expected_sample_ids = "['5a_HOM-I13', '5a_HOM-I14', '5a_HOM-I20', '5a_HOM-I21', '5a_HOM-I7', '1b_HOM-G58']"
    expected_filenames = "['HOM_20ind_17SNPs_last_10_samples.vcf.gz', 'HOM_20ind_17SNPs_first_10_samples.vcf.gz']"

    command = f"query tsv '{query_string}' --bucket-name {bucket}"
    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0

    assert query_string in result.stdout
    assert expected_sample_ids in result.stdout
    assert expected_filenames in result.stdout


def test_bcftools_pipe_query(CONSTANTS):
    """Test running a bcftools pipe query using the CLI."""
    bucket = CONSTANTS["QUERY_BUCKET"]
    tsv_filter = "Area:West of Ireland,Northern Portugal;"
    arg_command = "view -s SAMPLES; view -r 21:15000000-25000000"

    command = f"query bcftools-pipe --tsv-filter '{tsv_filter}' --command '{arg_command}' --bucket-name {bucket}"
    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0
    assert "Job submitted" in result.stdout

    time.sleep(3)  # Allow enough time for the job to complete
    command = f"files list --bucket-name {bucket}"
    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0
    assert "merged.vcf.gz" in result.stdout


@pytest.mark.parametrize(
    "bucket_name,tsv_filter,command,expected_error",
    [
        # Invalid bucket name
        ("non_existent_bucket", "Area:West of Ireland,Northern Portugal;", "DEFAULT", "ClientError"),
        # Malformed tsv filter (missing colon)
        ("DEFAULT", "Area West of Ireland", "DEFAULT", "SidecarInvalidFilterError"),
        # bad command
        ("DEFAULT", "DEFAULT", "invalid-command", "Unsupported bcftools command"),
        # empty command string
        ("DEFAULT", "DEFAULT", "", "Empty"),
    ],
)
def test_bcftools_pipe_query_errors(bucket_name, tsv_filter, command, expected_error, CONSTANTS):
    """
    Test bad formatted input raises errors

    TODO - these sorts of errors in the future should be handled by the API, and not sent to Celery.
    This test will need to be rewritten then, hence why it is quite sloppy now.
    """
    if "DEFAULT" in bucket_name:
        bucket_name = CONSTANTS["QUERY_BUCKET"]
    if "DEFAULT" in tsv_filter:
        tsv_filter = "Area:West of Ireland,Northern Portugal;"
    if "DEFAULT" in command:
        command = "view -s SAMPLES"

    command = f"query bcftools-pipe --tsv-filter '{tsv_filter}' --command '{command}' --bucket-name {bucket_name}"
    result = runner.invoke(app, shlex.split(command))

    task_id = result.stdout.strip().split()[-1]
    # TODO, assertion below should become 1, when API has the role of validating input.
    assert result.exit_code == 0
    time.sleep(1)

    job_status = runner.invoke(app, shlex.split(f"query task-status {task_id}"))
    assert job_status.exit_code == 0
    assert expected_error in job_status.stdout


def test_get_task_status_by_task_id(CONSTANTS):
    """Get the status of a task by its ID."""
    bucket = CONSTANTS["QUERY_BUCKET"]
    tsv_filter = "Area:West of Ireland,Northern Portugal;"
    arg_command = "view -s SAMPLES; view -r 21:15000000-25000000"

    command = f"query bcftools-pipe --tsv-filter '{tsv_filter}' --command '{arg_command}' --bucket-name {bucket}"
    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0
    task_id = result.stdout.strip().split()[-1]

    other_result = runner.invoke(app, shlex.split(command))
    assert other_result.exit_code == 0
    other_task_id = other_result.stdout.strip().split()[-1]

    command = f"query task-status {task_id}"
    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0
    assert "Task ID" in result.stdout
    assert task_id in result.stdout
    assert other_task_id not in result.stdout
