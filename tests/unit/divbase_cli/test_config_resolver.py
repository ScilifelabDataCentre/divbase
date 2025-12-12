from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from divbase_cli.cli_exceptions import AuthenticationError, ProjectNameNotSpecifiedError
from divbase_cli.config_resolver import (
    ensure_logged_in,
    resolve_divbase_api_url,
    resolve_download_dir,
    resolve_project,
)
from divbase_cli.user_config import ProjectConfig


def test_ensure_logged_in_success():
    """Test that ensure_logged_in returns the logged_in_url if the user is logged in."""
    mock_config = MagicMock()
    mock_config.logged_in_url = "https://example.com"

    with patch("divbase_cli.config_resolver.load_user_config", return_value=mock_config):
        result = ensure_logged_in(config_path=Path("/mock/path"))
        assert result == "https://example.com"


def test_ensure_logged_in_when_not_logged_in():
    """Test that ensure_logged_in raises AuthenticationError if the user is not logged in."""
    mock_config = MagicMock()
    mock_config.logged_in_url = None

    with (
        patch("divbase_cli.config_resolver.load_user_config", return_value=mock_config),
        pytest.raises(AuthenticationError, match="You are not logged in"),
    ):
        ensure_logged_in(config_path=Path("/mock/path"))


def test_ensure_logged_in_wrong_url():
    """Test that ensure_logged_in raises AuthenticationError if the logged_in_url does not match the desired_url."""
    mock_config = MagicMock()
    mock_config.logged_in_url = "https://example.com"
    with (
        patch("divbase_cli.config_resolver.load_user_config", return_value=mock_config),
        pytest.raises(AuthenticationError, match="You are not logged in to the correct DivBase URL"),
    ):
        ensure_logged_in(config_path=Path("/mock/path"), desired_url="https://wrong-url.com")


def test_resolve_project_with_explicit_project_name():
    """Test that resolve_project returns the correct project when the name is explicitly provided."""
    mock_project = ProjectConfig(name="test_project", divbase_url="https://example.com")
    mock_config = MagicMock()
    mock_config.project_info.return_value = mock_project

    with patch("divbase_cli.config_resolver.load_user_config", return_value=mock_config):
        result = resolve_project(project_name="test_project", config_path=Path("/mock/path"))
        assert result == mock_project
        mock_config.project_info.assert_called_once_with("test_project")


def test_resolve_project_with_default_project():
    """Test that resolve_project falls back to the default project if no name is provided."""
    mock_project = ProjectConfig(name="default_project", divbase_url="https://example.com")
    mock_config = MagicMock()
    mock_config.default_project = "default_project"
    mock_config.project_info.return_value = mock_project

    with patch("divbase_cli.config_resolver.load_user_config", return_value=mock_config):
        result = resolve_project(project_name=None, config_path=Path("/mock/path"))
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
        resolve_project(project_name=None, config_path=Path("/mock/path"))


def test_resolve_divbase_api_url_with_explicit_url():
    """Test that resolve_divbase_api_url returns the provided URL if explicitly given."""
    result = resolve_divbase_api_url(url="https://example.com", config_path=Path("/mock/path"))
    assert result == "https://example.com"


def test_resolve_divbase_api_url_with_no_url_given():
    """Test that resolve_divbase_api_url falls back to the default project's API URL."""
    mock_project = ProjectConfig(name="default_project", divbase_url="https://example.com")
    mock_config = MagicMock()
    mock_config.default_project = "default_project"
    mock_config.project_info.return_value = mock_project

    with patch("divbase_cli.config_resolver.load_user_config", return_value=mock_config):
        result = resolve_divbase_api_url(url=None, config_path=Path("/mock/path"))
        assert result == "https://example.com"
        mock_config.project_info.assert_called_once_with(name="default_project")


def test_resolve_divbase_api_url_no_url_and_no_default_project():
    """Test that resolve_divbase_api_url raises ValueError if no url given and no default project"""
    mock_config = MagicMock()
    mock_config.default_project = None

    with (
        patch("divbase_cli.config_resolver.load_user_config", return_value=mock_config),
        pytest.raises(ValueError, match="No default project is set in your user config"),
    ):
        resolve_divbase_api_url(url=None, config_path=Path("/mock/path"))


def test_resolve_download_dir_with_explicit_dir():
    """Test that resolve_download_dir returns the provided directory if explicitly given."""
    result = resolve_download_dir(download_dir="/mock/download", config_path=Path("/mock/path"))
    assert result == Path("/mock/download")


def test_resolve_download_dir_with_default_dir():
    """Test that resolve_download_dir falls back to the default directory in the user config."""
    mock_config = MagicMock()
    mock_config.default_download_dir = "/mock/default"

    with patch("divbase_cli.config_resolver.load_user_config", return_value=mock_config):
        result = resolve_download_dir(download_dir=None, config_path=Path("/mock/path"))
        assert result == Path("/mock/default")


def test_resolve_download_dir_with_current_dir():
    """Test that resolve_download_dir defaults to the current working directory if no directory is specified."""
    result = resolve_download_dir(download_dir=".", config_path=Path("/mock/path"))
    assert result == Path.cwd()
