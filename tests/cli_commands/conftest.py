"""
Setup for pytest fixtures for e2e testing of the CLI commands.

Pytest fixtures set up (and tear down) a test environment with the full DivBase stack running locally,

The S3FileManager class is patched in all tests to use this test MinIO server,
it is autoused, so it does not need to be specified in each test.
"""

import logging
from pathlib import Path

import pytest
from typer.testing import CliRunner

from divbase_cli.divbase_cli import app
from tests.helpers.minio_setup import (
    MINIO_URL,
)

runner = CliRunner()

logger = logging.getLogger(__name__)


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
    result = runner.invoke(app, create_command)

    assert result.exit_code == 0
    assert tmp_path.exists(), "Config file was not created"

    return fresh_config_path


@pytest.fixture
def user_config_path(tmp_path, CONSTANTS):
    """
    Fixture to provide a path to an "existing" user configuration file with
    some existing projects and a default project set.
    """
    existing_config_path = tmp_path / ".divbase_config.yaml"
    create_command = f"config create --config-file {existing_config_path}"
    result = runner.invoke(app, create_command)
    assert result.exit_code == 0

    for project in CONSTANTS["PROJECT_CONTENTS"]:
        add_command = f"config add-project {project} --divbase-url http://localhost:8001/api --s3-url {MINIO_URL} --config {existing_config_path}"
        result = runner.invoke(app, add_command)
        assert result.exit_code == 0

    set_default_command = f"config set-default {CONSTANTS['DEFAULT_PROJECT']} --config {existing_config_path}"
    result = runner.invoke(app, set_default_command)
    assert result.exit_code == 0

    return existing_config_path


@pytest.fixture
def fixtures_dir():
    """Path to the fixtures directory."""
    return Path(__file__).parent.parent / "fixtures"
