"""
Tests for the "divbase-cli dimensions" subcommand

NOTE: All tests are run against a MinIO server on localhost from docker-compose.
NOTE: The clean dimensions fixture ensures that the dimensions file is removed before and after each test,
"""

from unittest.mock import patch

import boto3
import pytest
from typer.testing import CliRunner

from divbase_tools.s3_client import create_s3_file_manager
from divbase_tools.tasks import update_vcf_dimensions_task
from divbase_tools.vcf_dimension_indexing import DIMENSIONS_FILE_NAME, VCFDimensionIndexManager
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


@patch("divbase_tools.tasks.create_s3_file_manager")
def test_update_vcf_dimensions_task_with_test_minio(mock_create_s3_manager, CONSTANTS):
    """
    Unit test that directly calls the update_vcf_dimensions_task task to create and update dimensions for the split-scaffold-project.
    Patches the S3 URL to use the test MinIO url.

    After asserting that the task logic has completed, create a new VCFDimensionIndexManager instance to read the .vcf_dimensions.yaml file
    with manager.get_dimensions_info() and assert that all VCF files in bucket have been indexed for their dimensions.
    """

    test_minio_url = CONSTANTS["MINIO_URL"]
    bucket_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]

    mock_create_s3_manager.side_effect = lambda url=None: create_s3_file_manager(url=test_minio_url)
    result = update_vcf_dimensions_task(bucket_name=bucket_name)
    assert result["status"] == "completed"

    s3_file_manager = create_s3_file_manager(url=test_minio_url)
    manager = VCFDimensionIndexManager(bucket_name=bucket_name, s3_file_manager=s3_file_manager)
    dimensions_info = manager.get_dimensions_info()
    vcf_files = [f for f in PROJECTS[bucket_name] if f.endswith(".vcf.gz") or f.endswith(".vcf")]
    indexed_files = [entry["filename"] for entry in dimensions_info.get("dimensions", [])]

    for vcf_file in vcf_files:
        assert vcf_file in indexed_files
