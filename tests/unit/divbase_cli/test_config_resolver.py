from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from pydantic import SecretStr

from divbase_cli.cli_config import cli_settings
from divbase_cli.cli_exceptions import AuthenticationError, ProjectNameNotSpecifiedError, ProjectNotInConfigError
from divbase_cli.config_resolver import (
    resolve_and_authenticate_project,
    resolve_download_dir,
    resolve_url_for_non_project_specific_commands,
)
from divbase_cli.user_config import ProjectConfig, UserConfig
from divbase_lib.divbase_constants import PAT_TOKEN_PREFIX

PROJECT_NAME = "test_project"
PROJECT_URL = "https://divbase.com"


@pytest.fixture(autouse=True)
def no_pat(monkeypatch):
    """
    Ensures DIVBASE_API_PAT is never taken from test runner environment (aka dev's token),
    this can interfere here.
    """
    monkeypatch.setattr(cli_settings, "DIVBASE_API_PAT", None)


@pytest.fixture
def mock_config_and_project():
    """
    Patches load_user_config to return a UserConfig with one project ("test_project")
    Tests can modify the mock_config, mock_project as needed
    """
    mock_project = ProjectConfig(name=PROJECT_NAME, divbase_url=PROJECT_URL)
    mock_config = UserConfig(config_path=Path("/mock/config.yaml"), projects=[mock_project])

    with patch("divbase_cli.config_resolver.load_user_config", return_value=mock_config):
        yield mock_config, mock_project


def test_resolve_and_authenticate_project_success(mock_config_and_project):
    """Test that resolve_and_authenticate_project returns the project config if the user is logged in to the right URL."""
    mock_config, mock_project = mock_config_and_project
    mock_config.logged_in_url = PROJECT_URL

    result = resolve_and_authenticate_project(project_name=PROJECT_NAME)
    assert result == mock_project


def test_resolve_and_authenticate_project_with_default_project(mock_config_and_project):
    """Test that resolve_and_authenticate_project falls back to the default project if no name is provided."""
    mock_config, mock_project = mock_config_and_project
    mock_config.default_project = PROJECT_NAME
    mock_config.logged_in_url = PROJECT_URL

    result = resolve_and_authenticate_project(project_name=None)
    assert result == mock_project


def test_resolve_and_authenticate_project_for_non_default_project(mock_config_and_project):
    """Test that resolve_and_authenticate_project can resolve a non default project if name is provided."""
    mock_config, mock_project = mock_config_and_project
    mock_config.default_project = PROJECT_NAME
    mock_config.logged_in_url = PROJECT_URL

    second_project = ProjectConfig(name="second_project", divbase_url=PROJECT_URL)
    mock_config.projects.append(second_project)
    result = resolve_and_authenticate_project(project_name="second_project")
    assert result == second_project


def test_resolve_and_authenticate_project_no_project_specified(mock_config_and_project):
    """Test that resolve_and_authenticate_project raises ProjectNameNotSpecifiedError if no project can be resolved."""
    mock_config, _ = mock_config_and_project
    mock_config.default_project = None

    with pytest.raises(ProjectNameNotSpecifiedError):
        resolve_and_authenticate_project(project_name=None)


def test_resolve_and_authenticate_project_not_in_config(mock_config_and_project):
    """Test that resolve_and_authenticate_project raises ProjectNotInConfigError if the named project isn't in the user's config."""
    mock_config, _ = mock_config_and_project
    unknown_project_name = "some_other_project"
    assert unknown_project_name not in mock_config.all_project_names

    with pytest.raises(ProjectNotInConfigError):
        resolve_and_authenticate_project(project_name=unknown_project_name)


def test_resolve_and_authenticate_project_wrong_url(mock_config_and_project):
    """Test that resolve_and_authenticate_project raises AuthenticationError if logged in to a different URL than the project's."""
    mock_config, _ = mock_config_and_project
    mock_config.logged_in_url = "https://wrong-url.com"

    with pytest.raises(AuthenticationError, match="You are trying to run a command"):
        resolve_and_authenticate_project(project_name=PROJECT_NAME)


def test_resolve_and_authenticate_project_with_pat(mock_config_and_project, monkeypatch):
    """resolve_and_authenticate_project returns the project config when no active session but a PAT is set."""
    mock_config, mock_project = mock_config_and_project
    mock_config.logged_in_url = None
    monkeypatch.setattr(cli_settings, "DIVBASE_API_PAT", SecretStr(f"{PAT_TOKEN_PREFIX}_abc123"))

    result = resolve_and_authenticate_project(project_name=PROJECT_NAME)
    assert result == mock_project


def test_resolve_and_authenticate_project_not_logged_in_and_no_pat(mock_config_and_project):
    """Test that resolve_and_authenticate_project raises AuthenticationError if not logged in and no PAT is set."""
    mock_config, _ = mock_config_and_project
    mock_config.logged_in_url = None

    with pytest.raises(AuthenticationError, match="You are not logged in"):
        resolve_and_authenticate_project(project_name=PROJECT_NAME)


def test_resolve_url_for_non_project_specific_commands_with_pat(monkeypatch):
    """resolve_url_for_non_project_specific_commands returns DIVBASE_API_URL when no active session and PAT is set."""
    mock_config = MagicMock()
    mock_config.logged_in_url = None
    monkeypatch.setattr(cli_settings, "DIVBASE_API_PAT", SecretStr(f"{PAT_TOKEN_PREFIX}_abc123"))
    monkeypatch.setattr(cli_settings, "DIVBASE_API_URL", "https://default.example.com")

    with patch("divbase_cli.config_resolver.load_user_config", return_value=mock_config):
        result = resolve_url_for_non_project_specific_commands()

    assert result == "https://default.example.com"


def test_resolve_url_for_non_project_specific_commands_without_pat_and_not_logged_in():
    """resolve_url_for_non_project_specific_commands raises AuthenticationError when not logged in and no PAT."""
    mock_config = MagicMock()
    mock_config.logged_in_url = None

    with (
        patch("divbase_cli.config_resolver.cli_settings") as mock_settings,
        patch("divbase_cli.config_resolver.load_user_config", return_value=mock_config),
        pytest.raises(AuthenticationError, match="You are not logged in"),
    ):
        mock_settings.DIVBASE_API_PAT = None
        resolve_url_for_non_project_specific_commands()


def test_resolve_url_for_non_project_specific_commands_without_pat_and_logged_in():
    """resolve_url_for_non_project_specific_commands returns logged_in_url when no PAT but user is logged in."""
    mock_config = MagicMock()
    mock_config.logged_in_url = "https://logged-in.example.com"

    with (
        patch("divbase_cli.config_resolver.cli_settings") as mock_settings,
        patch("divbase_cli.config_resolver.load_user_config", return_value=mock_config),
    ):
        mock_settings.DIVBASE_API_PAT = None
        result = resolve_url_for_non_project_specific_commands()

    assert result == "https://logged-in.example.com"


def test_resolve_download_dir_with_explicit_dir():
    """Test that resolve_download_dir returns the provided directory if explicitly given."""
    result = resolve_download_dir(
        download_dir="/mock/download",
    )
    assert result == Path("/mock/download")


def test_resolve_download_dir_with_default_dir():
    """Test that resolve_download_dir falls back to the default directory in the user config."""
    mock_config = MagicMock()
    mock_config.default_download_dir = "/mock/default"

    with patch("divbase_cli.config_resolver.load_user_config", return_value=mock_config):
        result = resolve_download_dir(
            download_dir=None,
        )
        assert result == Path("/mock/default")


def test_resolve_download_dir_with_current_dir():
    """Test that resolve_download_dir defaults to the current working directory if no directory is specified."""
    result = resolve_download_dir(
        download_dir=".",
    )
    assert result == Path.cwd()
