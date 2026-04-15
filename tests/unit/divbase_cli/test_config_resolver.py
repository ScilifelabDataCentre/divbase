from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from divbase_cli.cli_exceptions import AuthenticationError, ProjectNameNotSpecifiedError
from divbase_cli.config_resolver import (
    ensure_logged_in,
    resolve_download_dir,
    resolve_project,
    resolve_url_for_non_project_specific_commands,
)
from divbase_cli.user_config import ProjectConfig


def test_ensure_logged_in_success():
    """Test that ensure_logged_in returns the logged_in_url if the user is logged in."""
    mock_config = MagicMock()
    mock_config.logged_in_url = "https://example.com"

    with patch("divbase_cli.config_resolver.load_user_config", return_value=mock_config):
        result = ensure_logged_in()
        assert result == "https://example.com"


def test_ensure_logged_in_when_not_logged_in():
    """Test that ensure_logged_in raises AuthenticationError if the user is not logged in."""
    mock_config = MagicMock()
    mock_config.logged_in_url = None

    with (
        patch("divbase_cli.config_resolver.load_user_config", return_value=mock_config),
        pytest.raises(AuthenticationError, match="You are not logged in"),
    ):
        ensure_logged_in()


def test_ensure_logged_in_wrong_url():
    """Test that ensure_logged_in raises AuthenticationError if the logged_in_url does not match the desired_url."""
    mock_config = MagicMock()
    mock_config.logged_in_url = "https://example.com"
    with (
        patch("divbase_cli.config_resolver.load_user_config", return_value=mock_config),
        pytest.raises(AuthenticationError, match="You are not logged in to the correct DivBase URL"),
    ):
        ensure_logged_in(desired_url="https://wrong-url.com")


def test_ensure_logged_in_with_pat_and_desired_url():
    """ensure_logged_in returns desired_url immediately when DIVBASE_API_PAT is set, skipping the login check."""
    with patch("divbase_cli.config_resolver.cli_settings") as mock_settings:
        mock_settings.DIVBASE_API_PAT = "divbase_pat_abc123"
        mock_settings.DIVBASE_API_URL = "https://default.example.com"

        result = ensure_logged_in(desired_url="https://project.example.com")

    assert result == "https://project.example.com"


def test_ensure_logged_in_with_pat_and_no_desired_url():
    """ensure_logged_in returns DIVBASE_API_URL when PAT is set but no desired_url is given."""
    with patch("divbase_cli.config_resolver.cli_settings") as mock_settings:
        mock_settings.DIVBASE_API_PAT = "divbase_pat_abc123"
        mock_settings.DIVBASE_API_URL = "https://default.example.com"

        result = ensure_logged_in()

    assert result == "https://default.example.com"


def test_resolve_project_with_explicit_project_name():
    """Test that resolve_project returns the correct project when the name is explicitly provided."""
    mock_project = ProjectConfig(name="test_project", divbase_url="https://example.com")
    mock_config = MagicMock()
    mock_config.project_info.return_value = mock_project

    with patch("divbase_cli.config_resolver.load_user_config", return_value=mock_config):
        result = resolve_project(
            project_name="test_project",
        )
        assert result == mock_project
        mock_config.project_info.assert_called_once_with("test_project")


def test_resolve_project_with_default_project():
    """Test that resolve_project falls back to the default project if no name is provided."""
    mock_project = ProjectConfig(name="default_project", divbase_url="https://example.com")
    mock_config = MagicMock()
    mock_config.default_project = "default_project"
    mock_config.project_info.return_value = mock_project

    with patch("divbase_cli.config_resolver.load_user_config", return_value=mock_config):
        result = resolve_project(
            project_name=None,
        )
        assert result == mock_project
        mock_config.project_info.assert_called_once_with("default_project")


def test_resolve_project_no_project_specified():
    """Test that resolve_project raises ProjectNameNotSpecifiedError if no project is specified."""
    mock_config = MagicMock()
    mock_config.default_project = None

    with (
        patch("divbase_cli.config_resolver.load_user_config", return_value=mock_config),
        pytest.raises(ProjectNameNotSpecifiedError),
    ):
        resolve_project(
            project_name=None,
        )


def test_resolve_url_for_non_project_specific_commands_with_pat():
    """resolve_url_for_non_project_specific_commands returns DIVBASE_API_URL when PAT is set."""
    with patch("divbase_cli.config_resolver.cli_settings") as mock_settings:
        mock_settings.DIVBASE_API_PAT = "divbase_pat_abc123"
        mock_settings.DIVBASE_API_URL = "https://default.example.com"

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
