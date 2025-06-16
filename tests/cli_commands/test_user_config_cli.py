"""
Tests for the user configuration CLI commands.

Tests start from 1 of 3 types of user config files:
- fresh_config: A newly created config file with no buckets.
- existing_config: A config file with some pre-existing buckets and a default bucket set.
- No config file: To test the "config create" command...
"""

import shlex

import pytest
from typer.testing import CliRunner

from divbase_tools.divbase_cli import app
from divbase_tools.user_config import load_user_config

runner = CliRunner()

BUCKETS = ["bucket1", "bucket2", "bucket3"]
DEFAULT_BUCKET = "bucket1"


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
def existing_config(tmp_path):
    """
    Fixture to provide a path to an "existing" configuration file with:
    some existing buckets and a default bucket set.
    """
    existing_config_path = tmp_path / ".divbase_config.yaml"
    create_command = f"config create --config-file {existing_config_path}"
    result = runner.invoke(app, shlex.split(create_command))
    assert result.exit_code == 0

    for bucket in BUCKETS:
        add_command = f"config add-bucket {bucket} --config {existing_config_path}"
        result = runner.invoke(app, shlex.split(add_command))
        assert result.exit_code == 0
    runner.invoke(app, shlex.split(f"config set-default {DEFAULT_BUCKET} --config {existing_config_path}"))

    return existing_config_path


def test_create_config_command(tmp_config_path):
    command = f"config create --config-file {tmp_config_path}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0
    assert tmp_config_path.exists(), "Config file was not created at the temporary path"


def test_cant_create_config_if_exists(tmp_config_path):
    command = f"config create --config-file {tmp_config_path}"
    result1 = runner.invoke(app, shlex.split(command))
    assert result1.exit_code == 0

    result2 = runner.invoke(app, shlex.split(command))
    assert result2.exit_code != 0
    assert isinstance(result2.exception, FileExistsError)


def test_add_bucket_command(fresh_config):
    command = f"config add-bucket test_bucket --config {fresh_config}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0
    config_contents = load_user_config(fresh_config)
    assert "test_bucket" in config_contents["buckets"], "Bucket was not added to config file"


def test_add_bucket_as_default_command(fresh_config):
    command = f"config add-bucket test_bucket --default --config {fresh_config}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0
    config_contents = load_user_config(fresh_config)
    assert "test_bucket" in config_contents["buckets"]
    assert config_contents["default_bucket"] == "test_bucket"


def test_set_default_bucket_command(existing_config):
    config_contents = load_user_config(existing_config)
    assert config_contents["default_bucket"] != "bucket2"

    command = f"config set-default bucket2 --config {existing_config}"
    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0

    config_contents = load_user_config(existing_config)
    assert config_contents["default_bucket"] == "bucket2"


def test_set_default_a_bucket_that_does_not_exist(existing_config):
    """Test setting a default bucket that does not exist in the config."""
    config_contents = load_user_config(existing_config)
    assert config_contents["default_bucket"] != "bucket2"

    command = f"config set-default bucket2 --config {existing_config}"
    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0

    config_contents = load_user_config(existing_config)
    assert config_contents["default_bucket"] == "bucket2"


def test_show_default_bucket_command(existing_config):
    command = f"config show-default --config {existing_config}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0
    assert "bucket" in result.output


def test_show_default_with_no_default_set_command(fresh_config):
    command = f"config show-default --config {fresh_config}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0
    assert "No default bucket" in result.output


def test_show_user_config_command(existing_config):
    command = f"config show --config {existing_config}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0
    for bucket in BUCKETS:
        assert bucket in result.output
    assert f"default_bucket: {DEFAULT_BUCKET}" in result.output


def test_show_user_config_with_no_buckets_command(fresh_config):
    command = f"config show --config {fresh_config}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0
    assert "buckets: []" in result.output
    assert "default_bucket: None" in result.output


remove_commands = [
    (BUCKETS[0], True),
    (BUCKETS[1], True),
    (BUCKETS[2], True),
    ("bucket-that-does-not-exist", False),
]


@pytest.mark.parametrize("bucket_name,outcome", remove_commands)
def test_remove_bucket_command(existing_config, bucket_name, outcome):
    command = f"config remove-bucket {bucket_name} --config {existing_config}"
    result = runner.invoke(app, shlex.split(command))

    if outcome:
        assert result.exit_code == 0
        config_contents = load_user_config(existing_config)
        assert bucket_name not in config_contents["buckets"]

    # should fail if bucket does not exist
    else:
        assert result.exit_code != 0
        assert isinstance(result.exception, ValueError)
