"""
Tests for the user configuration CLI commands.

Tests start from 1 of 3 types of user config files:
- fresh_config: A newly created config file with no buckets.
- user_config_path: A config file with some pre-existing buckets and a default bucket set.
- No config file: To test the "config create" command...
"""

import shlex

import pytest
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
    bucket_name = "test_bucket"
    command = f"config add-bucket test_bucket --config {fresh_config}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0
    config_contents = load_user_config(fresh_config)
    assert bucket_name in config_contents.all_bucket_names


def test_add_bucket_as_default_command(fresh_config):
    bucket_name = "test_bucket"

    command = f"config add-bucket {bucket_name} --default --config {fresh_config}"
    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0

    config = load_user_config(fresh_config)
    assert bucket_name in config.all_bucket_names
    assert config.default_bucket == bucket_name


def test_add_bucket_and_specify_urls(fresh_config):
    bucket_name = "test_bucket"
    divbase_url = "https://divbasewebsite.com"
    s3_url = "http://s3.divbasewebsite.com"

    command = f"config add-bucket {bucket_name} --divbase-url {divbase_url} --s3-url {s3_url} --config {fresh_config}"
    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0

    config_contents = load_user_config(fresh_config)
    bucket = config_contents.bucket_info(name=bucket_name)

    assert bucket is not None
    assert bucket.name == bucket_name
    assert bucket.divbase_url == divbase_url
    assert bucket.s3_url == s3_url


def test_add_bucket_that_already_exists(fresh_config):
    """
    Should warn user and not be duplicated in the config.
    The new version of bucket should not be added.
    """
    bucket_name = "test_bucket"
    initial_divbase_url = "https://divbasewebsite.se"
    new_divbase_url = "https://newdivbasewebsite.se"

    initial_command = f"config add-bucket {bucket_name} --divbase-url {initial_divbase_url} --config {fresh_config}"
    result = runner.invoke(app, shlex.split(initial_command))
    assert result.exit_code == 0

    user_config = load_user_config(fresh_config)
    bucket_info = user_config.bucket_info(name=bucket_name)
    assert bucket_info.divbase_url == initial_divbase_url
    assert bucket_info.name in user_config.all_bucket_names

    new_command = f"config add-bucket {bucket_name} --divbase-url {new_divbase_url} --config {fresh_config}"

    with pytest.warns(UserWarning, match=f"The bucket: '{bucket_name}' already existed"):
        result = runner.invoke(app, shlex.split(new_command))
    assert result.exit_code == 0

    user_config = load_user_config(fresh_config)
    bucket_info = user_config.bucket_info(name=bucket_name)
    assert bucket_info.divbase_url == new_divbase_url
    assert bucket_info.name in user_config.all_bucket_names


def test_set_default_bucket_command(user_config_path):
    user_config = load_user_config(user_config_path)
    assert user_config.default_bucket != "bucket2"

    command = f"config set-default bucket2 --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0

    user_config = load_user_config(user_config_path)
    assert user_config.default_bucket == "bucket2"


def test_set_default_a_bucket_that_does_not_exist(user_config_path):
    """Test setting a default bucket that does not exist in the config."""
    command = f"config set-default made-up-bucket --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))
    assert isinstance(result.exception, ValueError)


def test_show_default_bucket_command(user_config_path, CONSTANTS):
    command = f"config show-default --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0
    assert CONSTANTS["DEFAULT_BUCKET"] in result.output


def test_show_default_with_no_default_set_command(fresh_config):
    command = f"config show-default --config {fresh_config}"
    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0
    assert "No default bucket" in result.output


def test_remove_bucket_command(user_config_path, CONSTANTS):
    for bucket_name in CONSTANTS["BUCKET_CONTENTS"]:
        command = f"config remove-bucket {bucket_name} --config {user_config_path}"
        result = runner.invoke(app, shlex.split(command))
        assert result.exit_code == 0

        user_config = load_user_config(user_config_path)
        assert bucket_name not in user_config.all_bucket_names


def test_remove_bucket_that_does_not_exist(user_config_path):
    bucket_name = "does-not-exist"
    command = f"config remove-bucket {bucket_name} --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0
    assert f"The bucket '{bucket_name}' was not found in your config file" in result.output


def test_remove_default_bucket_command(user_config_path, CONSTANTS):
    command = f"config remove-bucket {CONSTANTS['DEFAULT_BUCKET']} --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0

    user_config = load_user_config(user_config_path)
    assert CONSTANTS["DEFAULT_BUCKET"] not in user_config.all_bucket_names
    assert user_config.default_bucket is None


def test_set_default_dload_dir_command(user_config_path):
    download_dir = "/tmp/downloads"

    command = f"config set-dload-dir {download_dir} --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0

    user_config = load_user_config(user_config_path)
    assert user_config.default_download_dir == download_dir


def test_show_user_config_command(user_config_path, CONSTANTS):
    command = f"config show --config {user_config_path}"
    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0

    for bucket_name in CONSTANTS["BUCKET_CONTENTS"]:
        assert bucket_name in result.output


def test_show_user_config_with_no_buckets_command(fresh_config):
    command = f"config show --config {fresh_config}"
    result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0
    assert "No buckets" in result.output
