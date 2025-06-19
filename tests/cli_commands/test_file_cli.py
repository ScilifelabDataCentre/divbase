"""
Tests for the "divbase-cli files" commands

All tests are run against a MinIO server on localhost from docker-compose.
"""

import shlex

from typer.testing import CliRunner

from divbase_tools.divbase_cli import app
from divbase_tools.exceptions import FilesAlreadyInBucketError

runner = CliRunner()


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


def test_upload_with_safe_mode(user_config_path, CONSTANTS):
    # TODO - think about implementation,
    # ideally have a "clean" bucket or known not upload file (by any other test) to test this on?
    pass


def test_upload_with_safe_mode_prevents_upload(user_config_path, CONSTANTS, fixtures_dir):
    test_files = [(fixtures_dir / file_name).resolve() for file_name in CONSTANTS["FILES_TO_UPLOAD_DOWNLOAD"]]

    command = f"files upload {' '.join(map(str, test_files))} --config {user_config_path} --safe-mode"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code != 0
    assert isinstance(result.exception, FilesAlreadyInBucketError)
