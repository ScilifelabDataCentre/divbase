"""
Tests for the "divbase-cli version" subcommand

NOTE: All tests are run against a MinIO server on localhost from docker-compose.
NOTE: The clean versions fixture ensures that the versioning file is removed before and after each test,
"""

import shlex

import boto3
import pytest
from minio_setup import MINIO_FAKE_ACCESS_KEY, MINIO_FAKE_SECRET_KEY, URL
from typer.testing import CliRunner

from divbase_tools.bucket_versioning import VERSION_FILE_NAME
from divbase_tools.divbase_cli import app
from divbase_tools.exceptions import BucketVersioningFileAlreadyExistsError, BucketVersionNotFoundError

runner = CliRunner()

VERSION_1_NAME = "v1.0.0"
VERSION_2_NAME = "v2.0.0"
VERSION_3_NAME = "v3.0.0"


@pytest.fixture(autouse=True)
def clean_versions(user_config_path, CONSTANTS):
    """
    Remove the versioning file and create a new one before each test.
    Used in all tests in this module.
    """
    for bucket_name in CONSTANTS["BUCKET_CONTENTS"]:
        s3_client = boto3.client(
            "s3",
            endpoint_url=URL,
            aws_access_key_id=MINIO_FAKE_ACCESS_KEY,
            aws_secret_access_key=MINIO_FAKE_SECRET_KEY,
        )
        s3_client.delete_object(Bucket=bucket_name, Key=VERSION_FILE_NAME)

    command = f"version create --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0
    assert "Bucket versioning file created" in result.stdout

    yield


def test_create_version_file_fails_if_already_exists(user_config_path):
    command = f"version create --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code != 0
    assert isinstance(result.exception, BucketVersioningFileAlreadyExistsError)


def test_create_version_for_non_default_bucket(user_config_path, CONSTANTS):
    bucket_name = CONSTANTS["NON_DEFAULT_BUCKET"]
    command = f"version create --bucket-name {bucket_name} --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0
    assert f"Bucket versioning file created in bucket: '{bucket_name}'" in result.stdout


def test_add_version(user_config_path):
    command = f"version add {VERSION_1_NAME} --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0
    assert f"New version: '{VERSION_1_NAME}'" in result.stdout


def test_add_version_with_description(user_config_path):
    description = "Initial release"
    command = f'version add {VERSION_1_NAME} --description "{description}" --config {user_config_path}'
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0

    list_cmd = f"version list --config {user_config_path}"
    list_result = runner.invoke(app, shlex.split(list_cmd))

    assert VERSION_1_NAME in list_result.stdout
    assert description in list_result.stdout


def test_add_multiple_versions(user_config_path):
    versions = [VERSION_1_NAME, VERSION_2_NAME, VERSION_3_NAME]
    for version in versions:
        runner.invoke(app, shlex.split(f"version add {version} --config {user_config_path}"))

    command = f"version list --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0

    for version in versions:
        assert version in result.stdout


def test_list_versions(user_config_path):
    runner.invoke(app, shlex.split(f"version add {VERSION_1_NAME} --config {user_config_path}"))
    runner.invoke(app, shlex.split(f"version add {VERSION_2_NAME} --config {user_config_path}"))

    command = f"version list --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0
    assert f"- '{VERSION_1_NAME}':" in result.stdout
    assert f"- '{VERSION_2_NAME}':" in result.stdout


def test_delete_version(user_config_path):
    runner.invoke(app, shlex.split(f"version add {VERSION_1_NAME} --config {user_config_path}"))

    command = f"version delete {VERSION_1_NAME} --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0
    assert f"version: '{VERSION_1_NAME}' was deleted" in result.stdout

    list_cmd = f"version list --config {user_config_path}"
    list_result = runner.invoke(app, shlex.split(list_cmd))

    assert VERSION_1_NAME not in list_result.stdout


def test_delete_nonexistent_version(user_config_path):
    nonexistent_version = "v99.99.99"
    command = f"version delete {nonexistent_version} --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code != 0
    assert isinstance(result.exception, BucketVersionNotFoundError)


def test_get_version_info(user_config_path, CONSTANTS):
    runner.invoke(app, shlex.split(f"version add {VERSION_1_NAME} --config {user_config_path}"))

    command = f"version info {VERSION_1_NAME} --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0
    assert f"'{VERSION_1_NAME}'" in result.stdout
    assert f"'{CONSTANTS['DEFAULT_BUCKET']}'" in result.stdout

    CONSTANTS["BUCKET_CONTENTS"]

    print(result.stdout)


def test_get_version_updates_hashes_on_new_upload(user_config_path, CONSTANTS, fixtures_dir):
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
    upload_cmd = f"files upload {test_file_path} --config {user_config_path}"

    result = runner.invoke(app, shlex.split(f"version add {VERSION_3_NAME} --config {user_config_path}"))
    assert result.exit_code == 0

    result = runner.invoke(app, shlex.split(upload_cmd))
    assert result.exit_code == 0

    result = runner.invoke(app, shlex.split(f"version add {VERSION_1_NAME} --config {user_config_path}"))
    assert result.exit_code == 0

    result = runner.invoke(app, shlex.split(upload_cmd))
    assert result.exit_code == 0

    result = runner.invoke(app, shlex.split(f"version add {VERSION_2_NAME} --config {user_config_path}"))
    assert result.exit_code == 0

    info_result_1 = runner.invoke(app, shlex.split(f"version info {VERSION_1_NAME} --config {user_config_path}"))
    info_result_2 = runner.invoke(app, shlex.split(f"version info {VERSION_2_NAME} --config {user_config_path}"))

    hash_v1, hash_v2 = "", ""
    for line in info_result_1.stdout.splitlines():
        if test_file_name in line:
            hash_v1 = line.split(f"- '{test_file_name}' : ")[1].strip()
            break
    for line in info_result_2.stdout.splitlines():
        if test_file_name in line:
            hash_v2 = line.split(f"- '{test_file_name}' : ")[1].strip()
            break

    assert hash_v1 != "", "Hash for bucket version 1 not found"
    assert hash_v2 != "", "Hash for bucket version 2 not found"
    assert hash_v1 != hash_v2, "Hashes for the same file in different versions should be different"
