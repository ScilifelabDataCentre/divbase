"""
Tests for the user configuration CLI commands.

Tests start from 1 of 3 types of user config files:
- fresh_config: A newly created user config file with no projects.
- user_config_path: A config file with some pre-existing projects and a default project set.
- No config file: To test the "config create" command...
"""

import pytest
from typer.testing import CliRunner

from divbase_cli.divbase_cli import app
from divbase_cli.user_config import load_user_config

runner = CliRunner()


def test_create_config_command(tmp_config_path):
    command = f"config create --config-file {tmp_config_path}"
    result = runner.invoke(app, command)

    assert result.exit_code == 0
    assert tmp_config_path.exists(), "Config file was not created at the temporary path"


def test_cant_create_config_if_exists(tmp_config_path):
    command = f"config create --config-file {tmp_config_path}"
    result1 = runner.invoke(app, command)
    assert result1.exit_code == 0

    result2 = runner.invoke(app, command)
    assert result2.exit_code != 0
    assert isinstance(result2.exception, FileExistsError)


def test_add_project_command(fresh_config_path):
    project_name = "test_project"
    command = f"config add-project {project_name}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0

    config_contents = load_user_config(fresh_config_path)
    assert project_name in config_contents.all_project_names


def test_add_project_as_default_command(fresh_config_path):
    project_name = "test_project"

    command = f"config add-project {project_name} --default"
    result = runner.invoke(app, command)
    assert result.exit_code == 0

    config = load_user_config(fresh_config_path)
    assert project_name in config.all_project_names
    assert config.default_project == project_name


def test_add_project_and_specify_urls(fresh_config_path):
    project_name = "test_project"
    divbase_url = "https://divbasewebsite.com"
    s3_url = "http://s3.divbasewebsite.com"

    command = f"config add-project {project_name} --divbase-url {divbase_url} --s3-url {s3_url}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0

    config_contents = load_user_config(fresh_config_path)
    project = config_contents.project_info(name=project_name)

    assert project is not None
    assert project.name == project_name
    assert project.divbase_url == divbase_url
    assert project.s3_url == s3_url


def test_add_project_that_already_exists(fresh_config_path):
    """
    Should warn user and not be duplicated in the config.
    The newer project settings should be in the config file (i.e. overwrite the old ones).
    """
    project_name = "test_project"
    initial_divbase_url = "https://divbasewebsite.se"
    new_divbase_url = "https://newdivbasewebsite.se"

    initial_command = f"config add-project {project_name} --divbase-url {initial_divbase_url}"
    result = runner.invoke(app, initial_command)
    assert result.exit_code == 0

    user_config = load_user_config(fresh_config_path)
    project_info = user_config.project_info(name=project_name)
    assert project_info.divbase_url == initial_divbase_url
    assert project_info.name in user_config.all_project_names

    new_command = f"config add-project {project_name} --divbase-url {new_divbase_url}"

    with pytest.warns(UserWarning, match=f"The project: '{project_name}' already existed"):
        result = runner.invoke(app, new_command)
    assert result.exit_code == 0

    user_config = load_user_config(fresh_config_path)
    project_info = user_config.project_info(name=project_name)
    assert project_info.divbase_url == new_divbase_url
    assert project_info.name in user_config.all_project_names


def test_set_default_project_command(user_config_path):
    user_config = load_user_config(user_config_path)
    assert user_config.default_project != "project2"

    command = "config set-default project2"
    result = runner.invoke(app, command)
    assert result.exit_code == 0

    user_config = load_user_config(user_config_path)
    assert user_config.default_project == "project2"


def test_set_default_project_for_a_project_that_does_not_exist(user_config_path):
    """Test setting a default project that does not exist in the config."""
    command = "config set-default made-up-project"
    result = runner.invoke(app, command)
    assert isinstance(result.exception, ValueError)


def test_show_default_project_command(user_config_path, CONSTANTS):
    command = "config show-default"
    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert CONSTANTS["DEFAULT_PROJECT"] in result.output


def test_show_default_with_no_default_set_command(fresh_config_path):
    command = "config show-default"
    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert "No default project" in result.output


def test_remove_project_command(user_config_path, CONSTANTS):
    for project_name in CONSTANTS["PROJECT_CONTENTS"]:
        command = f"config remove-project {project_name}"
        result = runner.invoke(app, command)
        assert result.exit_code == 0

        user_config = load_user_config(user_config_path)
        assert project_name not in user_config.all_project_names


def test_remove_project_that_does_not_exist(user_config_path):
    project_name = "does-not-exist"
    command = f"config remove-project {project_name}"
    result = runner.invoke(app, command)

    assert result.exit_code == 0
    assert f"The project '{project_name}' was not found in your config file" in result.output


def test_remove_default_project_command(user_config_path, CONSTANTS):
    command = f"config remove-project {CONSTANTS['DEFAULT_PROJECT']}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0

    user_config = load_user_config(user_config_path)
    assert CONSTANTS["DEFAULT_PROJECT"] not in user_config.all_project_names
    assert user_config.default_project is None


def test_set_default_dload_dir_command(user_config_path):
    download_dir = "/tmp/downloads"

    command = f"config set-dload-dir {download_dir}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0

    user_config = load_user_config(user_config_path)
    assert user_config.default_download_dir == download_dir


def test_show_user_config_command(user_config_path, CONSTANTS):
    command = "config show"
    result = runner.invoke(app, command)
    assert result.exit_code == 0

    for project_name in CONSTANTS["PROJECT_CONTENTS"]:
        assert project_name in result.output


def test_show_user_config_with_no_projects_command(fresh_config_path):
    command = "config show"
    result = runner.invoke(app, command)

    assert result.exit_code == 0
    assert "No projects" in result.output
