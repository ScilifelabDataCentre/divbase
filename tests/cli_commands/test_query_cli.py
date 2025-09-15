"""
e2e tests for the "divbase-cli query" commands

All tests are run against a docker compose setup with the entire DivBase stack running locally.

A project (CONSTANTS["QUERY_PROJECT"]) is made available with input files for the tests.
"""

import time

import boto3
import pytest
from typer.testing import CliRunner

from divbase_tools.divbase_cli import app
from divbase_tools.exceptions import ProjectNotInConfigError

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
