"""
Tests for the user configuration CLI commands.

Tests start from 1 of 3 types of user config files:
- fresh_config: A newly created config file with no buckets.
- user_config_path: A config file with some pre-existing buckets and a default bucket set.
- No config file: To test the "config create" command...
"""

import shlex

from typer.testing import CliRunner

from divbase_tools.divbase_cli import app
from divbase_tools.user_config import load_user_config

runner = CliRunner()


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


def test_add_bucket_that_already_exists(fresh_config):
    """Don't want to raise error, but should not be duplicated."""
    command = f"config add-bucket test_bucket --config {fresh_config}"

    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0

    config_contents = load_user_config(fresh_config)
    assert "test_bucket" in config_contents["buckets"]

    result = runner.invoke(app, shlex.split(command))


def test_set_default_bucket_command(user_config_path):
    config_contents = load_user_config(user_config_path)
    assert config_contents["default_bucket"] != "bucket2"

    command = f"config set-default bucket2 --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0

    config_contents = load_user_config(user_config_path)
    assert config_contents["default_bucket"] == "bucket2"


def test_set_default_a_bucket_that_does_not_exist(user_config_path):
    """Test setting a default bucket that does not exist in the config."""
    config_contents = load_user_config(user_config_path)
    assert config_contents["default_bucket"] != "bucket2"

    command = f"config set-default bucket2 --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0

    config_contents = load_user_config(user_config_path)
    assert config_contents["default_bucket"] == "bucket2"


def test_show_default_bucket_command(user_config_path):
    command = f"config show-default --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0
    assert "bucket" in result.output


def test_show_default_with_no_default_set_command(fresh_config):
    command = f"config show-default --config {fresh_config}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0
    assert "No default bucket" in result.output


def test_show_user_config_command(user_config_path, CONSTANTS):
    command = f"config show --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0
    for bucket in CONSTANTS["BUCKET_CONTENTS"]:
        assert bucket in result.output
    assert f"default_bucket: {CONSTANTS['DEFAULT_BUCKET']}" in result.output


def test_show_user_config_with_no_buckets_command(fresh_config):
    command = f"config show --config {fresh_config}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0
    assert "buckets: []" in result.output
    assert "default_bucket: None" in result.output


def test_remove_bucket_command(user_config_path, CONSTANTS):
    for bucket_name in CONSTANTS["BUCKET_CONTENTS"]:
        command = f"config remove-bucket {bucket_name} --config {user_config_path}"
        result = runner.invoke(app, shlex.split(command))
        assert result.exit_code == 0
        config_contents = load_user_config(user_config_path)
        assert bucket_name not in config_contents["buckets"]


def test_remove_bucket_that_does_not_exist(user_config_path):
    command = f"config remove-bucket bucket-that-does-not-exist --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code != 0
    assert isinstance(result.exception, ValueError)


def test_remove_default_bucket_command(user_config_path, CONSTANTS):
    command = f"config remove-bucket {CONSTANTS['DEFAULT_BUCKET']} --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0
    config_contents = load_user_config(user_config_path)
    assert CONSTANTS["DEFAULT_BUCKET"] not in config_contents["buckets"]
    assert config_contents["default_bucket"] is None
