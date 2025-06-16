"""
Tests for the "divbase-cli files" commands
"""

import shlex
from unittest.mock import MagicMock

import pytest
from typer.testing import CliRunner

from divbase_tools.divbase_cli import app

runner = CliRunner()


@pytest.fixture
def CONSTANTS():
    return {
        "BUCKETS": ["bucket1", "bucket2", "bucket3"],
        "DEFAULT_BUCKET": "bucket1",
        "FILES_IN_BUCKET": [
            "file1.vcf.gz",
            "file2.vcf.gz",
            "file3.vcf.gz",
        ],
        "FILES_TO_DOWNLOAD": "file1.vcf.gz",
    }


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


@pytest.fixture
def mock_s3_manager(CONSTANTS):
    """Mock S3FileManager for testing, mocks return values of the S3 managers methods."""
    mock_manager = MagicMock()

    mock_manager.list_files.return_value = CONSTANTS["FILES_IN_BUCKET"]
    mock_manager.download_file.return_value = CONSTANTS["FILES_TO_DOWNLOAD"]

    return mock_manager
