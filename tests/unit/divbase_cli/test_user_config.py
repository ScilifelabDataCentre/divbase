from typing import Generator

import pytest

from divbase_cli.cli_config import cli_settings
from divbase_cli.user_config import ProjectNotInConfigError, create_user_config, load_user_config


@pytest.fixture(autouse=True)
def clean_up_user_config() -> Generator[None, None, None]:
    """Ensure the user config and tokens file do not exist before each test and are removed after each test."""
    cli_settings.CONFIG_PATH.unlink(missing_ok=True)
    cli_settings.TOKENS_PATH.unlink(missing_ok=True)
    yield
    cli_settings.CONFIG_PATH.unlink(missing_ok=True)
    cli_settings.TOKENS_PATH.unlink(missing_ok=True)


def test_config_auto_made_if_no_config_file() -> None:
    """Test that the any command that needs the user config will build it if it doesn't exist yet."""
    assert cli_settings.CONFIG_PATH.exists() is False
    load_user_config()
    assert cli_settings.CONFIG_PATH.exists() is True


def test_add_project_creates_new_project() -> None:
    """Test that adding a project creates a new project in the config with the correct information."""
    create_user_config(config_path=cli_settings.CONFIG_PATH)
    user_config = load_user_config()
    user_config.add_project(name="project1", divbase_url="https://example1.com", is_default=True)

    assert len(user_config.projects) == 1
    assert user_config.projects[0].name == "project1"
    assert user_config.projects[0].divbase_url == "https://example1.com"
    assert user_config.default_project == "project1"


def test_add_project_replaces_existing_project() -> None:
    """Test that adding a project with the same name replaces the existing one."""

    create_user_config(config_path=cli_settings.CONFIG_PATH)
    user_config = load_user_config()
    user_config.add_project(name="project1", divbase_url="https://old.com", is_default=False)
    assert len(user_config.projects) == 1
    assert user_config.projects[0].divbase_url == "https://old.com"
    assert user_config.default_project is None

    with pytest.warns(UserWarning):
        user_config.add_project(name="project1", divbase_url="https://new.com", is_default=True)
    assert len(user_config.projects) == 1
    assert user_config.projects[0].divbase_url == "https://new.com"
    assert user_config.default_project == "project1"


def test_remove_project_on_existing_project() -> None:
    """Test that removing an existing project removes it from the config."""

    create_user_config(config_path=cli_settings.CONFIG_PATH)
    user_config = load_user_config()
    user_config.add_project(name="project1", divbase_url="https://example1.com", is_default=True)
    removed_project = user_config.remove_project(name="project1")

    assert removed_project == "project1"
    assert len(user_config.projects) == 0
    assert user_config.default_project is None


def test_remove_project_does_nothing_if_project_not_found() -> None:
    """Test that removing a non-existent project does nothing."""

    create_user_config(config_path=cli_settings.CONFIG_PATH)
    user_config = load_user_config()
    user_config.add_project(name="project1", divbase_url="https://example1.com", is_default=True)

    removed_project = user_config.remove_project(name="non_existent_project")

    assert removed_project is None
    assert len(user_config.projects) == 1


def test_set_default_project_updates_default() -> None:
    """Test that the default project is updated correctly."""

    create_user_config(config_path=cli_settings.CONFIG_PATH)
    user_config = load_user_config()
    user_config.add_project(name="project1", divbase_url="https://example1.com", is_default=False)
    user_config.add_project(name="project2", divbase_url="https://example2.com", is_default=False)

    default_project = user_config.set_default_project(name="project2")

    assert default_project == "project2"
    assert user_config.default_project == "project2"

    default_project = user_config.set_default_project(name="project1")

    assert default_project == "project1"
    assert user_config.default_project == "project1"


def test_set_default_project_raises_error_if_project_not_found() -> None:
    """Test that setting a default project raises an error if the project does not exist."""
    create_user_config(config_path=cli_settings.CONFIG_PATH)
    user_config = load_user_config()

    with pytest.raises(ValueError, match="was not found in your config file"):
        user_config.set_default_project(name="non_existent_project")


def test_set_default_download_dir_updates_directory() -> None:
    """Test that the default download directory is updated correctly."""
    create_user_config(config_path=cli_settings.CONFIG_PATH)
    user_config = load_user_config()

    new_dir = "/new/download/dir"
    updated_dir = user_config.set_default_download_dir(download_dir=new_dir)

    assert updated_dir == new_dir
    assert user_config.default_download_dir == new_dir


def test_project_info_returns_correct_project() -> None:
    """Test that project_info returns the correct project configuration."""
    create_user_config(config_path=cli_settings.CONFIG_PATH)
    user_config = load_user_config()

    user_config.add_project(name="project1", divbase_url="https://example1.com", is_default=False)

    project = user_config.project_info(name="project1")

    assert project.name == "project1"
    assert project.divbase_url == "https://example1.com"


def test_project_info_raises_error_if_project_not_found() -> None:
    """Test that project_info raises an error if the project does not exist."""
    create_user_config(config_path=cli_settings.CONFIG_PATH)
    user_config = load_user_config()

    with pytest.raises(ProjectNotInConfigError):
        user_config.project_info(name="non_existent_project")


def test_set_login_status_updates_login_info() -> None:
    """Test that the login status is updated correctly."""
    create_user_config(config_path=cli_settings.CONFIG_PATH)
    user_config = load_user_config()

    assert user_config.logged_in_url is None
    assert user_config.logged_in_email is None

    user_config.set_login_status(url="https://example.com", email="user@example.com")

    assert user_config.logged_in_url == "https://example.com"
    assert user_config.logged_in_email == "user@example.com"
