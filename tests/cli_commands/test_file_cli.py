"""
Tests for the "divbase-cli files" commands

All tests are run against a MinIO server on localhost from docker-compose.

A clean bucket (auto emptied before each test) is available to test that require a clean state.
"""

import shlex

import boto3
import pytest
from minio_setup import MINIO_FAKE_ACCESS_KEY, MINIO_FAKE_SECRET_KEY, URL
from typer.testing import CliRunner

from divbase_tools.divbase_cli import app
from divbase_tools.exceptions import FilesAlreadyInBucketError

runner = CliRunner()


@pytest.fixture(autouse=True)
def start_with_clean_bucket(CONSTANTS):
    """
    For tests that require a clean bucket, this fixture will
    ensure that the CLEANED_BUCKET bucket is empty before running the test.

    Caution:
    If you modify the approach make sure your implementation does not just add delete markers.
    The files need to be actually deleted.
    """
    s3_resource = boto3.resource(
        "s3", endpoint_url=URL, aws_access_key_id=MINIO_FAKE_ACCESS_KEY, aws_secret_access_key=MINIO_FAKE_SECRET_KEY
    )

    # pylance does not understand boto3 resource returns types, hence ignore below
    bucket = s3_resource.Bucket(CONSTANTS["CLEANED_BUCKET"])  # type: ignore
    bucket.object_versions.delete()

    yield


def test_list_files(user_config_path, CONSTANTS):
    """Test basic usage of files list command."""
    command = f"files list --config {user_config_path}"

    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0

    default_bucket_name = CONSTANTS["DEFAULT_BUCKET"]

    for file in CONSTANTS["BUCKET_CONTENTS"][default_bucket_name]:
        assert file in result.stdout, f"File {file} not found in the output of the list_files command"


def test_list_non_default_bucket(user_config_path, CONSTANTS):
    """Test list files for the non-default bucket."""
    non_default_bucket = CONSTANTS["NON_DEFAULT_BUCKET"]
    files_in_bucket = CONSTANTS["BUCKET_CONTENTS"][non_default_bucket]

    command = f"files list --config {user_config_path} --bucket-name {non_default_bucket}"
    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0
    for file in files_in_bucket:
        assert file in result.stdout


def test_list_files_empty_bucket(user_config_path, CONSTANTS):
    """Test list files for an empty bucket."""
    command = f"files list --config {user_config_path} --bucket-name {CONSTANTS['EMPTY_BUCKET']}"
    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0
    assert "No files found" in result.stdout


def test_upload_1_file(user_config_path, CONSTANTS, fixtures_dir):
    """Test upload 1 file to the bucket."""
    test_file = (fixtures_dir / CONSTANTS["FILES_TO_UPLOAD_DOWNLOAD"][0]).resolve()

    command = f"files upload {test_file} --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0
    assert f"{str(test_file)}" in result.stdout


def test_upload_1_file_to_non_default_bucket(user_config_path, CONSTANTS, fixtures_dir):
    test_file = (fixtures_dir / CONSTANTS["FILES_TO_UPLOAD_DOWNLOAD"][0]).resolve()

    command = f"files upload {test_file} --config {user_config_path} --bucket-name {CONSTANTS['NON_DEFAULT_BUCKET']}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0
    assert f"{str(test_file)}" in result.stdout


def test_upload_multiple_files_at_once(user_config_path, CONSTANTS, fixtures_dir):
    test_files = [(fixtures_dir / file_name).resolve() for file_name in CONSTANTS["FILES_TO_UPLOAD_DOWNLOAD"]]

    command = f"files upload {' '.join(map(str, test_files))} --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0
    for file in test_files:
        assert f"{str(file)}" in result.stdout


def test_upload_dir_contents(user_config_path, CONSTANTS, fixtures_dir):
    """Test upload all files in a directory."""
    files = [x for x in fixtures_dir.glob("*") if x.is_file()]  # does not get subdirs

    command = f"files upload --upload-dir {fixtures_dir.resolve()} --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0
    for file in files:
        assert f"{file.resolve()}" in result.stdout


def test_upload_with_safe_mode(user_config_path, CONSTANTS, fixtures_dir):
    """Test upload with safe mode works first time, but fails on subsequent attempts."""
    file_name = CONSTANTS["FILES_TO_UPLOAD_DOWNLOAD"][0]
    file_path = f"{fixtures_dir}/{file_name}"
    command = (
        f"files upload {file_path} --safe-mode --bucket-name {CONSTANTS['CLEANED_BUCKET']} --config {user_config_path} "
    )

    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0

    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code != 0
    assert isinstance(result.exception, FilesAlreadyInBucketError)


def test_no_file_uploaded_if_some_duplicated_with_safe_mode(user_config_path, CONSTANTS, fixtures_dir):
    """
    Test that no files are uploaded with safe mode on,
    if at least 1 of the files trying to be uploaded already exists in the bucket.
    """
    test_files = [(fixtures_dir / file_name).resolve() for file_name in CONSTANTS["FILES_TO_UPLOAD_DOWNLOAD"]]

    # upload just 1 of the files first
    command = f"files upload {test_files[0]} --safe-mode --bucket-name {CONSTANTS['CLEANED_BUCKET']} --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0

    # none should be uploaded as the first one already exists
    command = f"files upload {' '.join(map(str, test_files))} --safe-mode --bucket-name {CONSTANTS['CLEANED_BUCKET']} --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code != 0
    assert isinstance(result.exception, FilesAlreadyInBucketError)

    command = f"files list --bucket-name {CONSTANTS['CLEANED_BUCKET']} --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0
    assert test_files[0].name in result.stdout
    for file in test_files[1:]:
        assert file.name not in result.stdout, f"File {file.name} was uploaded when it shouldn't have been."
