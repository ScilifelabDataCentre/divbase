import time
from unittest.mock import MagicMock, patch

import pytest
from pydantic import SecretStr

from divbase_cli.user_auth import (
    TokenData,
    check_existing_session,
)


@pytest.fixture
def mock_logged_out_config():
    """Fixture to provide a mock user config"""
    config = MagicMock()
    config.logged_in_url = None
    config.logged_in_email = None
    config.set_login_status = MagicMock()
    return config


@pytest.fixture
def mock_logged_in_config():
    """Fixture to provide a mock user config"""
    config = MagicMock()
    config.logged_in_url = "https://divbase.com"
    config.logged_in_email = "test@divbase.com"
    config.set_login_status = MagicMock()
    return config


@pytest.fixture
def mock_token_data():
    """Fixture to provide a mock token data"""
    return TokenData(
        access_token=SecretStr("access"),
        refresh_token=SecretStr("refresh"),
        access_token_expires_at=9999999999,
        refresh_token_expires_at=9999999999,
    )


def test_is_access_token_expired_method(mock_token_data):
    """Test that is_access_token_expired correctly identifies expired tokens."""
    assert mock_token_data.is_access_token_expired() is False

    mock_token_data.access_token_expires_at = int(time.time()) - 10  # Expired 10 seconds ago
    assert mock_token_data.is_access_token_expired() is True


def test_is_refresh_token_expired_method(mock_token_data):
    """Test that is_refresh_token_expired correctly identifies expired tokens."""
    assert mock_token_data.is_refresh_token_expired() is False

    mock_token_data.refresh_token_expires_at = int(time.time()) - 10  # Expired 10 seconds ago
    assert mock_token_data.is_refresh_token_expired() is True


@patch("divbase_cli.user_auth.load_user_tokens")
def test_check_existing_session_not_logged_in(mock_logged_out_config):
    """Test that check_existing_session returns None if not logged in."""

    result = check_existing_session("https://divbase.com", mock_logged_out_config)
    assert result is None


@patch("divbase_cli.user_auth.load_user_tokens")
def test_check_existing_session_when_logged_in_to_different_divbase_url(mock_logged_in_config, mock_token_data):
    """Test that check_existing_session returns None if logged in to a different divbase URL."""

    result = check_existing_session("https://otherdivbaseinstance.com", mock_logged_in_config)
    assert result is None


@patch("divbase_cli.user_auth.load_user_tokens")
def test_check_existing_session_when_logged_in(mock_load_user_tokens, mock_logged_in_config, mock_token_data):
    """Test that check_existing_session returns refresh token expiry if logged in."""
    mock_load_user_tokens.return_value = mock_token_data

    result = check_existing_session("https://divbase.com", mock_logged_in_config)
    assert result == mock_token_data.refresh_token_expires_at
    mock_load_user_tokens.assert_called_once()
