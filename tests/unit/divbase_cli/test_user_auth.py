import contextlib
import json
import stat
import time
from unittest.mock import MagicMock, patch

import keyring
import pytest
from keyring.errors import KeyringError, NoKeyringError, PasswordDeleteError
from pydantic import SecretStr

from divbase_cli.cli_config import cli_settings
from divbase_cli.user_auth import (
    TokenData,
    _delete_stored_jwts,
    check_existing_session,
    load_user_tokens,
)

# cli config set to use "divbase-cli-test" (via pytest.ini) as the keyring service name,
# so tokens and PATs stored in tests don't interfere with real credentials in dev computer with service name "divbase-cli"
_TEST_KEYRING_SERVICE = "divbase-cli-test"


@pytest.fixture(autouse=True)
def clean_up_keyring_entries():
    """Clean up any test keyring entries added after each test"""
    yield
    with contextlib.suppress(NoKeyringError, PasswordDeleteError):
        keyring.delete_password(service_name=_TEST_KEYRING_SERVICE, username=cli_settings.KEYRING_USERNAME)


@pytest.fixture
def mock_logged_out_user_config():
    """Fixture to provide a mock user config"""
    config = MagicMock()
    config.logged_in_url = None
    config.logged_in_email = None
    config.set_login_status = MagicMock()
    return config


@pytest.fixture
def mock_logged_in_user_config():
    config = MagicMock()
    config.logged_in_url = "https://divbase.com"
    config.logged_in_email = "test@divbase.com"
    config.set_login_status = MagicMock()
    return config


@pytest.fixture
def mock_token_data():
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


def test_check_existing_session_not_logged_in(mock_logged_out_user_config):
    """Test that check_existing_session returns None if not logged in."""

    result = check_existing_session("https://divbase.com", mock_logged_out_user_config)
    assert result is None


@patch("divbase_cli.user_auth.load_user_tokens")
def test_check_existing_session_when_logged_in_to_different_divbase_url(mock_logged_in_user_config):
    """Test that check_existing_session returns None if logged in to a different divbase URL."""

    result = check_existing_session("https://otherdivbaseinstance.com", mock_logged_in_user_config)
    assert result is None


@patch("divbase_cli.user_auth.load_user_tokens")
def test_check_existing_session_when_logged_in(mock_load_user_tokens, mock_logged_in_user_config, mock_token_data):
    """Test that check_existing_session returns refresh token expiry if logged in."""
    mock_load_user_tokens.return_value = mock_token_data

    result = check_existing_session("https://divbase.com", mock_logged_in_user_config)
    assert result == mock_token_data.refresh_token_expires_at
    mock_load_user_tokens.assert_called_once()


@patch("divbase_cli.user_auth.keyring.set_password", side_effect=NoKeyringError)
def test_dump_tokens_falls_back_to_file_when_no_keyring(mock_set_password, mock_token_data, tmp_path):
    """
    If afftempt to dump tokens with keyring fails, tokens should be written to a file with 0600 permissions."""
    output = tmp_path / ".secrets"

    mock_token_data.dump_tokens(output_path=output)

    assert output.exists()
    file_mode = stat.S_IMODE(output.stat().st_mode)
    assert file_mode == 0o600


@patch("divbase_cli.user_auth.keyring.set_password", side_effect=KeyringError("backend error"))
def test_dump_tokens_falls_back_to_file_on_keyring_error(mock_set_password, mock_token_data, tmp_path):
    """dump_tokens writes a 0600 file when an unexpected keyring error occurs."""
    output = tmp_path / ".secrets"

    mock_token_data.dump_tokens(output_path=output)

    assert output.exists()
    file_mode = stat.S_IMODE(output.stat().st_mode)
    assert file_mode == 0o600


@patch("divbase_cli.user_auth.keyring.get_password")
def test_load_user_tokens_reads_from_keyring(mock_get_password, mock_token_data, tmp_path):
    """load_user_tokens returns tokens from keyring when present."""
    mock_get_password.return_value = json.dumps(
        {
            "access_token": "access",
            "refresh_token": "refresh",
            "access_token_expires_at": mock_token_data.access_token_expires_at,
            "refresh_token_expires_at": mock_token_data.refresh_token_expires_at,
        }
    )

    result = load_user_tokens(token_path=tmp_path / ".secrets")

    assert result is not None
    assert result.access_token.get_secret_value() == "access"
    assert result.refresh_token.get_secret_value() == "refresh"


@patch("divbase_cli.user_auth.keyring.get_password", side_effect=NoKeyringError)
def test_load_user_tokens_falls_back_to_file_when_no_keyring(mock_get_password, mock_token_data, tmp_path):
    """load_user_tokens reads from file when keyring is unavailable."""
    token_file = tmp_path / ".secrets"
    token_file.write_text(
        f"access_token: access\n"
        f"refresh_token: refresh\n"
        f"access_token_expires_at: {mock_token_data.access_token_expires_at}\n"
        f"refresh_token_expires_at: {mock_token_data.refresh_token_expires_at}\n"
    )

    result = load_user_tokens(token_path=token_file)

    assert result is not None
    assert result.access_token.get_secret_value() == "access"


@patch("divbase_cli.user_auth.keyring.get_password", return_value=None)
def test_load_user_tokens_returns_none_when_no_tokens(mock_get_password, tmp_path):
    """load_user_tokens returns None when neither keyring nor file has tokens."""
    result = load_user_tokens(token_path=tmp_path / ".secrets")
    assert result is None


@patch("divbase_cli.user_auth.keyring.delete_password", side_effect=PasswordDeleteError)
def test_delete_stored_tokens_tolerates_missing_keyring_entry(mock_delete_password, tmp_path):
    """_delete_stored_jwts does not raise if the keyring entry doesn't exist."""
    token_file = tmp_path / ".secrets"
    _delete_stored_jwts(token_file)
    _delete_stored_jwts(token_file)
