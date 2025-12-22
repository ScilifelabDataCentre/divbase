from pathlib import Path

import pytest

from divbase_cli.user_config import ProjectNotInConfigError, create_user_config, load_user_config


def test_show_user_config_fails_if_no_config_file(tmp_path: Path) -> None:
    """Test that the show command fails if no config file exists."""

    with pytest.raises(FileNotFoundError):
        load_user_config(config_path=tmp_path / "non_existent_config.yaml")


def test_add_project_creates_new_project(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yaml"
    create_user_config(config_path)
    user_config = load_user_config(config_path)

    user_config.add_project(name="project1", divbase_url="https://example1.com", is_default=True)

    assert len(user_config.projects) == 1
    assert user_config.projects[0].name == "project1"
    assert user_config.projects[0].divbase_url == "https://example1.com"
    assert user_config.default_project == "project1"


def test_add_project_replaces_existing_project(tmp_path: Path) -> None:
    """Test that adding a project with the same name replaces the existing one."""
    config_path = tmp_path / "config.yaml"
    create_user_config(config_path)
    user_config = load_user_config(config_path)

    user_config.add_project(name="project1", divbase_url="https://old.com", is_default=False)
    assert len(user_config.projects) == 1
    assert user_config.projects[0].divbase_url == "https://old.com"
    assert user_config.default_project is None

    with pytest.warns(UserWarning):
        user_config.add_project(name="project1", divbase_url="https://new.com", is_default=True)
    assert len(user_config.projects) == 1
    assert user_config.projects[0].divbase_url == "https://new.com"
    assert user_config.default_project == "project1"


def test_remove_project_on_existing_project(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yaml"
    create_user_config(config_path)
    user_config = load_user_config(config_path)

    user_config.add_project(name="project1", divbase_url="https://example1.com", is_default=True)
    removed_project = user_config.remove_project(name="project1")

    assert removed_project == "project1"
    assert len(user_config.projects) == 0
    assert user_config.default_project is None


def test_remove_project_does_nothing_if_project_not_found(tmp_path: Path) -> None:
    """Test that removing a non-existent project does nothing."""
    config_path = tmp_path / "config.yaml"
    create_user_config(config_path)
    user_config = load_user_config(config_path)

    user_config.add_project(name="project1", divbase_url="https://example1.com", is_default=True)

    removed_project = user_config.remove_project(name="non_existent_project")

    assert removed_project is None
    assert len(user_config.projects) == 1


def test_set_default_project_updates_default(tmp_path: Path) -> None:
    """Test that the default project is updated correctly."""
    config_path = tmp_path / "config.yaml"
    create_user_config(config_path)
    user_config = load_user_config(config_path)

    user_config.add_project(name="project1", divbase_url="https://example1.com", is_default=False)
    user_config.add_project(name="project2", divbase_url="https://example2.com", is_default=False)

    default_project = user_config.set_default_project(name="project2")

    assert default_project == "project2"
    assert user_config.default_project == "project2"

    default_project = user_config.set_default_project(name="project1")

    assert default_project == "project1"
    assert user_config.default_project == "project1"


def test_set_default_project_raises_error_if_project_not_found(tmp_path: Path) -> None:
    """Test that setting a default project raises an error if the project does not exist."""
    config_path = tmp_path / "config.yaml"
    create_user_config(config_path)
    user_config = load_user_config(config_path)

    with pytest.raises(ValueError, match="was not found in your config file"):
        user_config.set_default_project(name="non_existent_project")


def test_set_default_download_dir_updates_directory(tmp_path: Path) -> None:
    """Test that the default download directory is updated correctly."""
    config_path = tmp_path / "config.yaml"
    create_user_config(config_path)
    user_config = load_user_config(config_path)

    new_dir = "/new/download/dir"
    updated_dir = user_config.set_default_download_dir(download_dir=new_dir)

    assert updated_dir == new_dir
    assert user_config.default_download_dir == new_dir


def test_project_info_returns_correct_project(tmp_path: Path) -> None:
    """Test that project_info returns the correct project configuration."""
    config_path = tmp_path / "config.yaml"
    create_user_config(config_path)
    user_config = load_user_config(config_path)

    user_config.add_project(name="project1", divbase_url="https://example1.com", is_default=False)

    project = user_config.project_info(name="project1")

    assert project.name == "project1"
    assert project.divbase_url == "https://example1.com"


def test_project_info_raises_error_if_project_not_found(tmp_path: Path) -> None:
    """Test that project_info raises an error if the project does not exist."""
    config_path = tmp_path / "config.yaml"
    create_user_config(config_path)
    user_config = load_user_config(config_path)

    with pytest.raises(ProjectNotInConfigError):
        user_config.project_info(name="non_existent_project")


def test_set_login_status_updates_login_info(tmp_path: Path) -> None:
    """Test that the login status is updated correctly."""
    config_path = tmp_path / "config.yaml"
    create_user_config(config_path)
    user_config = load_user_config(config_path)

    assert user_config.logged_in_url is None
    assert user_config.logged_in_email is None

    user_config.set_login_status(url="https://example.com", email="user@example.com")

    assert user_config.logged_in_url == "https://example.com"
    assert user_config.logged_in_email == "user@example.com"
