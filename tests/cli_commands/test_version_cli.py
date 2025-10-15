"""
Tests for the "divbase-cli version" subcommand

NOTE: All tests are run against a MinIO server on localhost from docker-compose.
NOTE: The clean versions fixture ensures that the versioning file is removed before and after each test,
"""

import boto3
import pytest
from typer.testing import CliRunner

from divbase_cli.divbase_cli import app
from divbase_lib.exceptions import DivBaseAPIError

runner = CliRunner()


VERSION_1_NAME = "v1.0.0"
VERSION_2_NAME = "v2.0.0"
VERSION_3_NAME = "v3.0.0"

VERSION_FILE_NAME = ".bucket_versions.yaml"


@pytest.fixture(autouse=True)
def clean_versions(logged_in_edit_user_with_existing_config, CONSTANTS):
    """
    Remove the versioning file and create a new one before each test.
    Used in all tests in this module.
    """
    for project_name in CONSTANTS["PROJECT_CONTENTS"]:
        s3_client = boto3.client(
            "s3",
            endpoint_url=CONSTANTS["MINIO_URL"],
            aws_access_key_id=CONSTANTS["BAD_ACCESS_KEY"],
            aws_secret_access_key=CONSTANTS["BAD_SECRET_KEY"],
        )
        s3_client.delete_object(Bucket=project_name, Key=VERSION_FILE_NAME)

    command = "version create"
    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert "Bucket versioning file created" in result.stdout

    yield


def test_create_version_file_fails_if_already_exists(logged_in_edit_user_with_existing_config):
    command = "version create"
    result = runner.invoke(app, command)
    assert result.exit_code != 0
    assert isinstance(result.exception, DivBaseAPIError)
    assert result.exception.error_type == "bucket_versioning_file_already_exists_error"
    assert result.exception.status_code == 400


def test_create_version_for_non_default_project(logged_in_edit_user_with_existing_config, CONSTANTS):
    project_name = CONSTANTS["NON_DEFAULT_PROJECT"]
    command = f"version create --project {project_name}"
    result = runner.invoke(app, command)

    assert result.exit_code == 0
    assert f"Bucket versioning file created for project: '{project_name}'" in result.stdout


def test_add_version(logged_in_edit_user_with_existing_config):
    command = f"version add {VERSION_1_NAME}"
    result = runner.invoke(app, command)

    assert result.exit_code == 0
    assert f"New version: '{VERSION_1_NAME}'" in result.stdout


def test_add_version_with_description(logged_in_edit_user_with_existing_config):
    description = "Initial release"
    command = f'version add {VERSION_1_NAME} --description "{description}"'
    result = runner.invoke(app, command)

    assert result.exit_code == 0

    list_cmd = "version list"
    list_result = runner.invoke(app, list_cmd)

    assert VERSION_1_NAME in list_result.stdout
    assert description in list_result.stdout


def test_add_multiple_versions(logged_in_edit_user_with_existing_config):
    versions = [VERSION_1_NAME, VERSION_2_NAME, VERSION_3_NAME]
    for version in versions:
        command = f"version add {version}"
        runner.invoke(app, command)

    command = "version list"
    result = runner.invoke(app, command)

    assert result.exit_code == 0

    for version in versions:
        assert version in result.stdout


def test_attempt_add_version_that_already_exists_fails(logged_in_edit_user_with_existing_config):
    command = "version add v1.0.0"

    result = runner.invoke(app, command)
    assert result.exit_code == 0

    result = runner.invoke(app, command)
    assert result.exit_code != 0
    assert isinstance(result.exception, DivBaseAPIError)
    assert result.exception.error_type == "bucket_version_already_exists_error"
    assert result.exception.status_code == 400


def test_add_version_works_with_clean_project(logged_in_edit_user_with_existing_config, CONSTANTS):
    """
    Using the clean project which has no files at all in the bucket (not even the bucket metadata file)
    So validating you don't need to run `version create` first.
    """
    clean_project = CONSTANTS["CLEANED_PROJECT"]

    command = f"version add {VERSION_1_NAME} --project {clean_project}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0

    command = f"version list --project {clean_project}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert VERSION_1_NAME in result.stdout


def test_list_versions(logged_in_edit_user_with_existing_config):
    runner.invoke(app, f"version add {VERSION_1_NAME}")
    runner.invoke(app, f"version add {VERSION_2_NAME}")

    command = "version list"
    result = runner.invoke(app, command)

    assert result.exit_code == 0
    assert f"{VERSION_1_NAME}" in result.stdout
    assert f"{VERSION_2_NAME}" in result.stdout


def test_delete_version(logged_in_edit_user_with_existing_config):
    command = f"version add {VERSION_1_NAME}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0

    command = f"version delete {VERSION_1_NAME}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert f"version: '{VERSION_1_NAME}' was deleted" in result.stdout

    list_cmd = "version list"
    result = runner.invoke(app, list_cmd)
    assert result.exit_code == 0
    assert VERSION_1_NAME not in result.stdout


def test_delete_nonexistent_version(logged_in_edit_user_with_existing_config):
    nonexistent_version = "v99.99.99"
    command = f"version delete {nonexistent_version}"
    result = runner.invoke(app, command)

    assert result.exit_code != 0
    assert isinstance(result.exception, DivBaseAPIError)
    assert result.exception.error_type == "bucket_version_not_found_error"
    assert result.exception.status_code == 404


def test_get_version_info(logged_in_edit_user_with_existing_config, CONSTANTS):
    default_project = CONSTANTS["DEFAULT_PROJECT"]
    files_in_project = CONSTANTS["PROJECT_CONTENTS"][default_project]

    command = f"version add {VERSION_1_NAME}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0

    command = f"version info {VERSION_1_NAME}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0

    assert f"{VERSION_1_NAME}" in result.stdout
    assert f"{default_project}" in result.stdout

    for file in files_in_project:
        assert f"- '{file}' :" in result.stdout


def test_get_version_info_for_version_that_does_not_exist(logged_in_edit_user_with_existing_config, CONSTANTS):
    command = "version info does_not_exist"
    result = runner.invoke(app, command)

    assert result.exit_code != 0
    assert isinstance(result.exception, DivBaseAPIError)
    assert result.exception.error_type == "bucket_version_not_found_error"
    assert result.exception.status_code == 404


def test_get_version_updates_hashes_on_new_upload(logged_in_edit_user_with_existing_config, CONSTANTS, fixtures_dir):
    """
    Test a simple protocol where:
    1. upload file
    2. create a version
    3. upload the same file again
    4. create another version

    Validate that the hashes for the uploaded file in each version are different.
    """
    test_file_name = CONSTANTS["FILES_TO_UPLOAD_DOWNLOAD"][0]
    test_file_path = (fixtures_dir / test_file_name).resolve()
    upload_cmd = f"files upload {test_file_path}"

    result = runner.invoke(app, f"version add {VERSION_3_NAME}")
    assert result.exit_code == 0

    result = runner.invoke(app, upload_cmd)
    assert result.exit_code == 0

    result = runner.invoke(app, f"version add {VERSION_1_NAME}")
    assert result.exit_code == 0

    result = runner.invoke(app, upload_cmd)
    assert result.exit_code == 0

    result = runner.invoke(app, f"version add {VERSION_2_NAME}")
    assert result.exit_code == 0

    info_result_1 = runner.invoke(app, f"version info {VERSION_1_NAME}")
    assert info_result_1.exit_code == 0
    info_result_2 = runner.invoke(app, f"version info {VERSION_2_NAME}")
    assert info_result_2.exit_code == 0

    hash_v1, hash_v2 = "", ""
    for line in info_result_1.stdout.splitlines():
        if test_file_name in line:
            hash_v1 = line.split(f"- '{test_file_name}' : ")[1].strip()
            break
    for line in info_result_2.stdout.splitlines():
        if test_file_name in line:
            hash_v2 = line.split(f"- '{test_file_name}' : ")[1].strip()
            break

    assert hash_v1 != "", "Hash for version 1 not found"
    assert hash_v2 != "", "Hash for version 2 not found"
    assert hash_v1 != hash_v2, "Hashes for the same file in different versions should be different"
