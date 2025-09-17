"""
Tests for the "divbase-cli dimensions" subcommand

NOTE: All tests are run against a MinIO server on localhost from docker-compose.
NOTE: The clean dimensions fixture ensures that the dimensions file is removed before and after each test,
"""

import ast
import re
from unittest.mock import patch

import boto3
import pytest
import yaml
from typer.testing import CliRunner

from divbase_cli.divbase_cli import app
from divbase_lib.exceptions import NoVCFFilesFoundError, VCFDimensionsFileMissingOrEmptyError
from divbase_lib.s3_client import create_s3_file_manager
from divbase_lib.vcf_dimension_indexing import DIMENSIONS_FILE_NAME, VCFDimensionIndexManager
from divbase_worker.tasks import update_vcf_dimensions_task
from tests.helpers.minio_setup import PROJECTS

runner = CliRunner()


@pytest.fixture(autouse=True)
def clean_dimensions(user_config_path, CONSTANTS):
    """
    Remove the dimensions file and create a new one before each test.
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

    yield


def test_update_vcf_dimensions_task_directly(CONSTANTS, run_update_dimensions):
    """
    Test that first runs the run_update_dimensions fixture to create a .vcf_dimensions.yaml file in the bucket.
    Then it creates a new VCFDimensionIndexManager instance to read the .vcf_dimensions.yaml file
    with manager.get_dimensions_info() and assert that all VCF files in bucket have been indexed for their dimensions.
    """
    test_minio_url = CONSTANTS["MINIO_URL"]
    bucket_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    run_update_dimensions(bucket_name=bucket_name)

    s3_file_manager = create_s3_file_manager(url=test_minio_url)
    manager = VCFDimensionIndexManager(bucket_name=bucket_name, s3_file_manager=s3_file_manager)
    dimensions_info = manager.get_dimensions_info()
    vcf_files = [f for f in PROJECTS[bucket_name] if f.endswith(".vcf.gz") or f.endswith(".vcf")]
    indexed_files = [entry["filename"] for entry in dimensions_info.get("dimensions", [])]

    for vcf_file in vcf_files:
        assert vcf_file in indexed_files


def test_show_vcf_dimensions_task(CONSTANTS, run_update_dimensions, user_config_path):
    """
    Test that first runs the run_update_dimensions fixture to create a .vcf_dimensions.yaml file in the bucket.
    Then it runs the CLI command to show the VCF dimensions and asserts that the output is as expected.
    """

    bucket_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    run_update_dimensions(bucket_name=bucket_name)

    # Basic version of command
    command = f"dimensions show --project {bucket_name} --config {user_config_path}"
    cli_result = runner.invoke(app, command)
    assert cli_result.exit_code == 0

    yaml.safe_load(cli_result.stdout)

    vcf_files = [f for f in PROJECTS[bucket_name] if f.endswith(".vcf.gz") or f.endswith(".vcf")]
    for vcf_file in vcf_files:
        assert vcf_file in cli_result.stdout, f"{vcf_file} not found in CLI output:\n{cli_result.stdout}"

    # Unique-scaffolds version of command
    command = f"dimensions show --project {bucket_name} --config {user_config_path} --unique-scaffolds"
    cli_result = runner.invoke(app, command)
    assert cli_result.exit_code == 0

    lines = cli_result.stdout.splitlines()
    scaffold_line = next((line for line in lines if line.startswith("['")), None)
    assert scaffold_line is not None, f"Scaffold list not found in output:\n{cli_result.stdout}"

    scaffold_names = ast.literal_eval(scaffold_line)
    expected_scaffolds = ["1", "4", "5", "6", "7", "8", "13", "18", "20", "21", "22", "24"]
    assert scaffold_names == expected_scaffolds, f"Expected {expected_scaffolds}, got {scaffold_names}"

    # Filename version of command
    for vcf_file in vcf_files:
        command = f"dimensions show --project {bucket_name} --config {user_config_path} --filename {vcf_file}"
        cli_result = runner.invoke(app, command)
        assert cli_result.exit_code == 0
        match = re.search(r"HOM_20ind_17SNPs\.(\d+)\.vcf\.gz", vcf_file)
        scaffold_name = match.group(1)
        entry = yaml.safe_load(cli_result.stdout)
        scaffolds = entry.get("dimensions", {}).get("scaffolds", [])
        assert scaffold_name in scaffolds, f"{scaffold_name} not found in scaffolds for {vcf_file}: {scaffolds}"


def test_show_vcf_dimensions_task_when_file_missing(CONSTANTS, user_config_path, caplog):
    """
    Test runs CLI command to show the VCF dimensions handles the case when the dimensions file is missing.
    """

    bucket_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]

    command = f"dimensions show --project {bucket_name} --config {user_config_path}"

    result = runner.invoke(app, command)
    assert result.exit_code != 0
    assert isinstance(result.exception, VCFDimensionsFileMissingOrEmptyError)
    assert bucket_name in str(result.exception)


def test_get_dimensions_info_returns_empty(CONSTANTS):
    """
    Test that asserts that get_dimensions_info returns {"dimensions": []} when the dimensions file is missing or empty.
    The clean_dimensions fixture is run before this test, thus there should be no .vcf_dimensions.yaml file present.
    """
    test_minio_url = CONSTANTS["MINIO_URL"]
    bucket_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]

    s3_file_manager = create_s3_file_manager(url=test_minio_url)
    manager = VCFDimensionIndexManager(bucket_name=bucket_name, s3_file_manager=s3_file_manager)

    result = manager.get_dimensions_info()
    assert result == {"dimensions": []}


def test_update_vcf_dimensions_task_raises_no_vcf_files_error(CONSTANTS):
    """
    Test that runs the update_vcf_dimensions_task with a bucket that has no VCF files.
    It should raise an error that there are no VCF files to process.
    """
    test_minio_url = CONSTANTS["MINIO_URL"]
    bucket_name = "empty-project"

    with patch("divbase_worker.tasks.create_s3_file_manager") as mock_create_s3_manager:
        mock_create_s3_manager.side_effect = lambda url=None: create_s3_file_manager(url=test_minio_url)
        with pytest.raises(NoVCFFilesFoundError):
            update_vcf_dimensions_task(bucket_name=bucket_name)


def test_remove_VCF_and_update_dimension_entry(CONSTANTS):
    """
    Test that mocks removing a VCF file from the bucket, and updates the dimension entry based on the removal.
    """
    s3_client = boto3.client(
        "s3",
        endpoint_url=CONSTANTS["MINIO_URL"],
        aws_access_key_id=CONSTANTS["BAD_ACCESS_KEY"],
        aws_secret_access_key=CONSTANTS["BAD_SECRET_KEY"],
    )
    bucket_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    vcf_file = "HOM_20ind_17SNPs.8.vcf.gz"

    with patch.object(s3_client, "delete_object", return_value=None) as mock_delete:
        s3_client.delete_object(Bucket=bucket_name, Key=vcf_file)
        mock_delete.assert_called_once_with(Bucket=bucket_name, Key=vcf_file)

    s3_file_manager = create_s3_file_manager(url=CONSTANTS["MINIO_URL"])
    manager = VCFDimensionIndexManager(bucket_name=bucket_name, s3_file_manager=s3_file_manager)
    manager.remove_dimension_entry(vcf_file)
    filenames = [entry["filename"] for entry in manager.dimensions_info.get("dimensions", [])]
    assert vcf_file not in filenames


def test_update_vcf_dimensions_task_upload_failed(CONSTANTS):
    """
    Test that mocks that update_vcf_dimensions_task fails to upload the dimensions file to bucket.
    """

    bucket_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]

    with patch("divbase_worker.tasks.create_s3_file_manager") as mock_create_s3_manager:
        # Ensure test compose stack is used when running the task
        mock_create_s3_manager.side_effect = lambda url=None: create_s3_file_manager(url="http://localhost:9002")

        with patch(
            "divbase_lib.s3_client.S3FileManager.upload_str_as_s3_object",
            side_effect=Exception("Simulated upload failure"),
        ):
            result = update_vcf_dimensions_task(bucket_name=bucket_name, user_name="Test User")
            assert result["status"] == "error"
            assert "Failed to upload bucket dimensions file" in result["error"]
