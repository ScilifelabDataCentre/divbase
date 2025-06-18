"""
Tests for the "divbase-cli files" commands

Pytest fixtures set up (and tear down) a MinIO server on localhost with test buckets and files.

The S3FileManager class is patched in all tests to use this test MinIO server,
it is autoused, so it does not need to be specified in each test.

"""

import shlex
from unittest.mock import patch

import pytest
from minio_setup import BUCKETS, setup_minio_data, start_minio, stop_minio
from typer.testing import CliRunner

from divbase_tools.divbase_cli import app
from divbase_tools.s3_client import S3FileManager

runner = CliRunner()


@pytest.fixture(scope="session")
def CONSTANTS():
    return {
        "BAD_ACCESS_KEY": "minioadmin",
        "BAD_SECRET_KEY": "badpassword",
        "BUCKETS": BUCKETS,
        "DEFAULT_BUCKET": "bucket1",
        "NON_DEFAULT_BUCKET": "bucket2",
        "EMPTY_BUCKET": "empty-bucket",
        "FILES_TO_UPLOAD_DOWNLOAD": ["file1.txt", "file2.txt", "file3.txt"],
    }


@pytest.fixture(scope="session")
def minio_server():
    """Start Minio server, set up test buckets, users etc... and stop after all tests run"""
    try:
        start_minio()
        setup_minio_data()
        yield "http://localhost:9000"
    finally:
        print("Stopping Minio server...")
        stop_minio()


@pytest.fixture(autouse=True)
def patch_s3_file_manager(minio_server, CONSTANTS):
    """Fixture to patch create_s3_file_manager to use the test Minio server."""

    def mock_create_s3_file_manager():
        return S3FileManager(
            url=minio_server,
            access_key=CONSTANTS["BAD_ACCESS_KEY"],
            secret_key=CONSTANTS["BAD_SECRET_KEY"],
        )

    with patch("divbase_tools.services.create_s3_file_manager", side_effect=mock_create_s3_file_manager) as mock:
        yield mock


@pytest.fixture
def tmp_config_path(tmp_path):
    """
    Fixture to provide a path to where the a configuration can be created.
    """
    return tmp_path / "test_config.yaml"


@pytest.fixture
def fresh_config(tmp_path):
    """
    Fixture to provide a path to a pre-existing configuration file.
    """
    fresh_config_path = tmp_path / ".divbase_config.yaml"
    create_command = f"config create --config-file {fresh_config_path}"
    result = runner.invoke(app, shlex.split(create_command))

    assert result.exit_code == 0
    assert tmp_path.exists(), "Config file was not created"

    return fresh_config_path


@pytest.fixture
def user_config_path(tmp_path, CONSTANTS):
    """
    Fixture to provide a path to an "existing" user configuration file with
    some existing buckets and a default bucket set.
    """
    existing_config_path = tmp_path / ".divbase_config.yaml"
    create_command = f"config create --config-file {existing_config_path}"
    result = runner.invoke(app, shlex.split(create_command))
    assert result.exit_code == 0

    for bucket in CONSTANTS["BUCKETS"]:
        add_command = f"config add-bucket {bucket} --config {existing_config_path}"
        result = runner.invoke(app, shlex.split(add_command))
        assert result.exit_code == 0
    runner.invoke(app, shlex.split(f"config set-default {CONSTANTS['DEFAULT_BUCKET']} --config {existing_config_path}"))

    return existing_config_path
