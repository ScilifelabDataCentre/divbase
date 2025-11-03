"""
Tests for the "divbase-cli dimensions" subcommand
"""

import ast
import gzip
import os
import re
from unittest.mock import patch

import httpx
import pytest
import yaml
from typer.testing import CliRunner

from divbase_api.worker.tasks import update_vcf_dimensions_task
from divbase_cli.divbase_cli import app
from divbase_lib.exceptions import DivBaseAPIError, NoVCFFilesFoundError, VCFDimensionsFileMissingOrEmptyError
from divbase_lib.s3_client import create_s3_file_manager
from divbase_lib.vcf_dimension_indexing import create_vcf_dimension_manager
from tests.helpers.minio_setup import PROJECTS

runner = CliRunner()

api_base_url = os.environ["DIVBASE_API_URL"]


@pytest.fixture(autouse=True)
def clean_dimensions(logged_out_user_with_existing_config, CONSTANTS):
    """
    Clean up VCF dimensions from the database before each test using VCFDimensionIndexManager.
    """
    login_response = httpx.post(
        f"{api_base_url}/v1/auth/login",
        data={
            "grant_type": "password",
            "username": os.environ["WORKER_SERVICE_EMAIL"],
            "password": os.environ["WORKER_SERVICE_PASSWORD"],
        },
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        timeout=10.0,
    )
    login_response.raise_for_status()
    auth_token = login_response.json()["access_token"]

    for project_name in CONSTANTS["PROJECT_CONTENTS"]:
        try:
            dimension_manager = create_vcf_dimension_manager(bucket_name=project_name, auth_token=auth_token)
            vcf_dimensions_data = dimension_manager.get_dimensions_info()

            if vcf_dimensions_data.get("vcf_files"):
                project_id = dimension_manager._project_id
                for entry in vcf_dimensions_data["vcf_files"]:
                    vcf_file = entry["vcf_file_s3_key"]
                    try:
                        dimension_manager.delete_vcf_metadata(vcf_file, project_id)
                    except Exception as e:
                        print(f"Warning: Failed to delete VCF metadata for {vcf_file}: {e}")

            try:
                skipped_files = dimension_manager.get_skipped_files()
                if skipped_files and dimension_manager._project_id:
                    project_id = dimension_manager._project_id
                    for vcf_file in skipped_files:
                        try:
                            dimension_manager.delete_skipped_vcf(vcf_file, project_id)
                        except Exception as e:
                            print(f"Warning: Failed to delete skipped VCF entry for {vcf_file}: {e}")
            except Exception as e:
                print(f"Warning: Failed to clean skipped files for {project_name}: {e}")

        except ValueError:
            pass
        except Exception as e:
            print(f"Error cleaning up dimensions for {project_name}: {e}")

    yield


def test_update_vcf_dimensions_task_directly(CONSTANTS, run_update_dimensions):
    """
    Test that runs the update task and verifies all VCF files are indexed via the API.
    """
    bucket_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    result = run_update_dimensions(bucket_name=bucket_name)

    vcf_files = [f for f in PROJECTS[bucket_name] if f.endswith(".vcf.gz") or f.endswith(".vcf")]
    indexed_files = result.get("VCF files that were added to dimensions index by this job", [])

    for vcf_file in vcf_files:
        assert vcf_file in indexed_files, f"{vcf_file} not found in indexed files: {indexed_files}"


def test_show_vcf_dimensions_task(CONSTANTS, run_update_dimensions, logged_in_edit_user_with_existing_config):
    """
    Test the CLI show command after indexing dimensions via the API.
    """
    bucket_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    run_update_dimensions(bucket_name=bucket_name)

    # Basic version of command
    command = f"dimensions show --project {bucket_name}"
    cli_result = runner.invoke(app, command)
    assert cli_result.exit_code == 0

    dimensions_info = yaml.safe_load(cli_result.stdout)
    assert isinstance(dimensions_info, dict), f"Expected dict, got: {type(dimensions_info)}"
    indexed_files = dimensions_info.get("indexed_files", [])

    vcf_files = [f for f in PROJECTS[bucket_name] if f.endswith(".vcf.gz") or f.endswith(".vcf")]
    found_files = [entry.get("filename") for entry in indexed_files]
    for vcf_file in vcf_files:
        assert vcf_file in found_files, f"{vcf_file} not found in CLI output:\n{cli_result.stdout}"

    # Unique-scaffolds version of command
    command = f"dimensions show --project {bucket_name} --unique-scaffolds"
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
        command = f"dimensions show --project {bucket_name} --filename {vcf_file}"
        cli_result = runner.invoke(app, command)
        assert cli_result.exit_code == 0
        entry = yaml.safe_load(cli_result.stdout)
        assert entry.get("filename") == vcf_file
        match = re.search(r"HOM_20ind_17SNPs\.(\d+)\.vcf\.gz", vcf_file)
        if match:
            scaffold_name = match.group(1)
            scaffolds = entry.get("dimensions", {}).get("scaffolds", [])
            assert scaffold_name in scaffolds, f"{scaffold_name} not found in scaffolds for {vcf_file}: {scaffolds}"


def test_show_vcf_dimensions_task_when_file_missing(CONSTANTS, logged_in_edit_user_with_existing_config):
    """
    Test that the CLI handles the case when no dimensions are indexed (empty database).
    """
    bucket_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]

    command = f"dimensions show --project {bucket_name}"

    result = runner.invoke(app, command)
    assert result.exit_code != 0
    assert isinstance(result.exception, (VCFDimensionsFileMissingOrEmptyError, DivBaseAPIError))
    assert bucket_name in str(result.exception)


def test_get_dimensions_info_returns_empty(CONSTANTS):
    """
    Test that get_dimensions_info returns empty when no dimensions are indexed in the database.
    """
    bucket_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]

    api_base_url = "http://localhost:8001/api"
    login_response = httpx.post(
        f"{api_base_url}/v1/auth/login",
        data={
            "grant_type": "password",
            "username": os.environ["WORKER_SERVICE_EMAIL"],
            "password": os.environ["WORKER_SERVICE_PASSWORD"],
        },
        headers={"Content-Type": "application/x-www-form-urlencoded"},
    )
    login_response.raise_for_status()
    auth_token = login_response.json()["access_token"]

    manager = create_vcf_dimension_manager(bucket_name=bucket_name, auth_token=auth_token)
    result = manager.get_dimensions_info()
    assert result == {"vcf_files": []}  # TODO this might need to be updated


def test_update_vcf_dimensions_task_raises_no_vcf_files_error(CONSTANTS):
    """
    Test that the update task raises an error when the bucket has no VCF files.
    """
    test_minio_url = CONSTANTS["MINIO_URL"]
    bucket_name = "empty-project"

    with patch("divbase_api.worker.tasks.create_s3_file_manager") as mock_create_s3_manager:
        mock_create_s3_manager.side_effect = lambda url=None: create_s3_file_manager(url=test_minio_url)
        with pytest.raises(NoVCFFilesFoundError):
            update_vcf_dimensions_task(bucket_name=bucket_name)


def test_remove_VCF_and_update_dimension_entry(CONSTANTS):
    """
    Test removing a VCF metadata entry via the manager.
    """
    bucket_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    vcf_file = "HOM_20ind_17SNPs.8.vcf.gz"

    api_base_url = "http://localhost:8000/api"
    login_response = httpx.post(
        f"{api_base_url}/v1/auth/login",
        data={
            "grant_type": "password",
            "username": os.environ["WORKER_SERVICE_EMAIL"],
            "password": os.environ["WORKER_SERVICE_PASSWORD"],
        },
        headers={"Content-Type": "application/x-www-form-urlencoded"},
    )
    login_response.raise_for_status()
    auth_token = login_response.json()["access_token"]

    manager = create_vcf_dimension_manager(bucket_name=bucket_name, auth_token=auth_token)

    if manager._project_id:
        manager.delete_vcf_metadata(vcf_file, manager._project_id)
        updated_dimensions = manager.get_dimensions_info()
        filenames = [entry["vcf_file_s3_key"] for entry in updated_dimensions.get("vcf_files", [])]
        assert vcf_file not in filenames


@patch("divbase_api.worker.tasks.create_s3_file_manager")
def test_update_dimensions_skips_divbase_generated_vcf(mock_create_s3_manager, CONSTANTS, tmp_path):
    """
    Test that after running a query (which generates a DivBase result VCF),
    update_vcf_dimensions_task skips that file and returns a skip message.
    """
    mock_create_s3_manager.side_effect = lambda url=None: create_s3_file_manager(url=CONSTANTS["MINIO_URL"])

    bucket_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]

    divbase_vcf_name = "merged_test_divbase_result.vcf.gz"
    vcf_path = tmp_path / divbase_vcf_name
    divbase_vcf_content = (
        "##fileformat=VCFv4.2\n"
        '##DivBase_created="This is a results file created by a DivBase query; Date=Mon Oct 20 12:00:00 2025"\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n"
        "1\t1000\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n"
    )

    with gzip.open(vcf_path, "wt") as f:
        f.write(divbase_vcf_content)

    s3_file_manager = create_s3_file_manager(url=CONSTANTS["MINIO_URL"])

    s3_file_manager.upload_files(
        to_upload={divbase_vcf_name: vcf_path},
        bucket_name=bucket_name,
    )

    result = update_vcf_dimensions_task(bucket_name=bucket_name, user_name="Test User")

    skipped_files = result.get("VCF files skipped by this job (previous DivBase-generated result VCFs)", [])
    assert any(divbase_vcf_name in msg for msg in skipped_files), (
        f"Expected that this file was skipped; {divbase_vcf_name}; got: {skipped_files}"
    )

    assert result["status"] == "completed"


@patch("divbase_api.worker.tasks.create_s3_file_manager")
def test_update_dimensions_twice_with_no_new_VCF_added_inbetween(mock_create_s3_manager, CONSTANTS, tmp_path):
    """
    Test that after running update_vcf_dimensions_task twice with no new VCF files added in between,
    the task returns a message indicating no new files were found.
    """
    mock_create_s3_manager.side_effect = lambda url=None: create_s3_file_manager(url=CONSTANTS["MINIO_URL"])

    bucket_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]

    result_first_run = update_vcf_dimensions_task(bucket_name=bucket_name, user_name="Test User")

    assert result_first_run["status"] == "completed"
    added_files = result_first_run["VCF files that were added to dimensions index by this job"]
    expected_files = PROJECTS["split-scaffold-project"]
    expected_vcfs = [f for f in expected_files if f.endswith(".vcf.gz")]
    for vcf in expected_vcfs:
        assert vcf in added_files, f"{vcf} not found in indexed files: {added_files}"

    result_second_run = update_vcf_dimensions_task(bucket_name=bucket_name, user_name="Test User")
    assert result_second_run["status"] == "completed"
    assert (
        "None: no new VCF files or file versions were detected in the project."
        in result_second_run["VCF files that were added to dimensions index by this job"]
    )
