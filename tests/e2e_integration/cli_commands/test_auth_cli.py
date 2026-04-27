"""
E2E tests for the "divbase-cli auth" CLI commands.
"""

import shutil

import keyring
import pytest
from keyring.errors import NoKeyringError
from pydantic import SecretStr
from typer.testing import CliRunner

from divbase_cli.cli_config import cli_settings
from divbase_cli.cli_exceptions import AuthenticationError, DivBaseAPIConnectionError, DivBaseAPIError
from divbase_cli.divbase_cli import app
from divbase_cli.user_auth import LOGIN_AGAIN_MESSAGE

runner = CliRunner()

USER_EMAIL = "edit@divbase.se"
USER_PASSWORD = "badpassword"


@pytest.fixture()
def disable_keyring_backend(monkeypatch):
    """
    Some tests want to manipulate the JWTs (to e.g. simulate making them being expired or revoked).
    Awkward to do that if stored in device keyring. So we can disable keyring for these tests,
    which will fall back to storing the JWTs in a file.

    In e.g. CI where a keyring backend won't be available, tests will still work as they will from the start use the file-based fallback.
    """
    monkeypatch.setattr(keyring, "set_password", lambda *a, **kw: (_ for _ in ()).throw(NoKeyringError()))
    monkeypatch.setattr(keyring, "get_password", lambda *a, **kw: None)


def make_tokens_expired(access: bool = False, refresh: bool = False):
    """
    Helper function to make either the access token or refresh token expired in the tokens file.
    Sets the expiry time to 1970 (unix time stamp)
    """
    with open(cli_settings.TOKENS_PATH, "r") as token_file:
        lines = token_file.readlines()
    with open(cli_settings.TOKENS_PATH, "w") as token_file:
        for line in lines:
            if access and line.startswith("access_token_expires_at:"):
                token_file.write("access_token_expires_at: 1\n")
            elif refresh and line.startswith("refresh_token_expires_at:"):
                token_file.write("refresh_token_expires_at: 1\n")
            else:
                token_file.write(line)


def log_in_as_user():
    """
    Helper function to log in as a user.
    Not a fixture as want to ensure logged out before and after each test and could be timing issues when
    combined with other fixtures (e.g. ensure_logged_out).
    """
    command = f"auth login {USER_EMAIL}"
    result = runner.invoke(app=app, args=command, input=f"{USER_PASSWORD}\n")
    assert result.exit_code == 0
    assert "Logged in successfully" in result.stdout
    assert USER_EMAIL in result.stdout


def test_login_command():
    log_in_as_user()


def test_login_command_fails_with_invalid_credentials():
    """Test login command fails with invalid credentials."""
    command = f"auth login {USER_EMAIL}"

    result = runner.invoke(app=app, args=command, input="wrongpassword\n")
    assert result.exit_code != 0
    assert isinstance(result.exception, AuthenticationError)
    assert "Invalid email or password" in str(result.exception)


def test_login_command_with_invalid_server_url():
    """Test login command fails with an invalid server URL."""
    command = f"auth login {USER_EMAIL} --divbase-url https://invalid-url"

    result = runner.invoke(app=app, args=command, input=f"{USER_PASSWORD}\n")
    assert result.exit_code != 0
    assert isinstance(result.exception, DivBaseAPIConnectionError)


def test_login_command_already_logged_in():
    """Test login command when already logged in."""
    log_in_as_user()

    command = f"auth login {USER_EMAIL}"

    # cancel login, when warned already logged in
    result = runner.invoke(app=app, args=command, input=f"{USER_PASSWORD}\nn\n")
    assert result.exit_code == 0
    assert "Already logged in to" in result.stdout
    assert "Do you want to login again?" in result.stdout
    assert "Login cancelled." in result.stdout

    # don't cancel login process when warned already logged in
    result = runner.invoke(app=app, args=command, input=f"{USER_PASSWORD}\ny\n")
    assert result.exit_code == 0
    assert "Already logged in to" in result.stdout
    assert "Do you want to login again?" in result.stdout
    assert "Logged in successfully" in result.stdout
    assert USER_EMAIL in result.stdout


def test_force_login_option():
    """Should not prompt about logging in again"""
    log_in_as_user()
    command = f"auth login {USER_EMAIL} --force"

    result = runner.invoke(app=app, args=command, input=f"{USER_PASSWORD}\n")
    assert result.exit_code == 0
    assert "Logged in successfully" in result.stdout
    assert USER_EMAIL in result.stdout


def test_logout_command():
    """Test basic usage of logout and that running multiple times is ok."""
    command = "auth logout"

    result = runner.invoke(app=app, args=command)
    assert result.exit_code == 0
    assert "Logged out successfully" in result.stdout

    result = runner.invoke(app=app, args=command)
    assert result.exit_code == 0
    assert "Logged out successfully" in result.stdout


def test_login_logout_cycle():
    """Test a few repeated login/logout cycles."""
    login_command = f"auth login {USER_EMAIL} --force"
    logout_command = "auth logout"

    for _ in range(3):
        result = runner.invoke(app=app, args=login_command, input=f"{USER_PASSWORD}\n")
        assert result.exit_code == 0
        assert "Logged in successfully" in result.stdout
        assert USER_EMAIL in result.stdout

        result = runner.invoke(app=app, args=logout_command)
        assert result.exit_code == 0
        assert "Logged out successfully" in result.stdout


def test_whoami_command():
    """Test basic usage of whoami command."""
    log_in_as_user()
    command = "auth whoami"

    result = runner.invoke(app=app, args=command)
    assert result.exit_code == 0
    assert USER_EMAIL in result.stdout


def test_whoami_command_fails_if_not_logged_in():
    """Test basic usage of whoami command."""
    command = "auth whoami"

    result = runner.invoke(app=app, args=command)
    assert result.exit_code != 0
    assert isinstance(result.exception, AuthenticationError)


def test_whoami_command_needing_refresh_token(disable_keyring_backend):
    """
    We simulate that the access token has expired by manually setting it to be expired in the users .secrets file

    This will force the CLI to use the refresh token to get a new access token when running the whoami command.
    """
    log_in_as_user()
    make_tokens_expired(access=True)

    whoami_command = "auth whoami"
    result = runner.invoke(app=app, args=whoami_command)
    assert result.exit_code == 0
    assert USER_EMAIL in result.stdout


def test_whoami_with_expired_tokens_fails(disable_keyring_backend):
    """
    Simulate that both the access token and refresh token have expired by manually setting them to be expired
    in the users .secrets file.

    This should force authentication to fail when running the whoami command, requiring the user to log in again.
    """
    log_in_as_user()
    make_tokens_expired(access=True, refresh=True)

    command = "auth whoami"
    result = runner.invoke(app=app, args=command)
    assert result.exit_code != 0
    assert isinstance(result.exception, AuthenticationError)
    assert LOGIN_AGAIN_MESSAGE in str(result.exception)


def test_using_revoked_refresh_token_fails(tmp_path, disable_keyring_backend):
    """
    Validate that when the refresh token is revoked, the user cannot use it to get a new access token.

    To do this:
    - Log in as user
    - copy tokens + config file
    - log out (this revokes the refresh token, server side
    - copy back in old tokens + config file
    - make access token expired
    - run whoami command, should fail because revoked refresh token.
    - log in again, should work
    - run whoami command, should work
    """
    logout_command = "auth logout"
    whoami_command = "auth whoami"

    log_in_as_user()

    shutil.copy(cli_settings.CONFIG_PATH, tmp_path / "config_backup.yaml")
    shutil.copy(cli_settings.TOKENS_PATH, tmp_path / "tokens_backup.yaml")

    result = runner.invoke(app=app, args=logout_command)
    assert result.exit_code == 0

    shutil.copy(tmp_path / "tokens_backup.yaml", cli_settings.TOKENS_PATH)
    shutil.copy(tmp_path / "config_backup.yaml", cli_settings.CONFIG_PATH)
    make_tokens_expired(access=True)

    result = runner.invoke(app=app, args=whoami_command)
    assert result.exit_code != 0
    assert isinstance(result.exception, AuthenticationError)
    assert LOGIN_AGAIN_MESSAGE in str(result.exception)

    log_in_as_user()

    result = runner.invoke(app=app, args=whoami_command)
    assert result.exit_code == 0
    assert USER_EMAIL in result.stdout


def test_login_with_outdated_cli_version_fails(monkeypatch):
    """Test that login fails if the CLI version is outdated (rejected by the API middleware)"""
    monkeypatch.setattr("divbase_cli.user_auth.cli_version", "0.0.0")

    command = f"auth login {USER_EMAIL} --force"
    result = runner.invoke(app=app, args=command, input=f"{USER_PASSWORD}\n")

    assert result.exit_code != 0
    assert isinstance(result.exception, DivBaseAPIError)
    assert "400" in str(result.exception)
    assert "cli_version_outdated_error" in str(result.exception)


def test_any_command_with_outdated_cli_version_fails(monkeypatch):
    """Test that any command fails if the CLI version is outdated."""
    log_in_as_user()

    monkeypatch.setattr("divbase_cli.user_auth.cli_version", "0.0.0")

    command = "auth whoami"
    result = runner.invoke(app=app, args=command)

    assert result.exit_code != 0
    assert isinstance(result.exception, DivBaseAPIError)
    assert "400" in str(result.exception)
    assert "cli_version_outdated_error" in str(result.exception)


def test_jwt_session_takes_priority_over_pat(disable_keyring_backend, monkeypatch):
    """
    When a JWT session (stored in fallback file) and a PAT are both available, the JWT is used.

    We know it works as if the PAT was used, would get error as PAT made up.
    """
    log_in_as_user()

    monkeypatch.setattr(cli_settings, "DIVBASE_API_PAT", SecretStr("divbase_pat_fake_not_a_real_token"))

    result = runner.invoke(app=app, args="auth whoami")
    assert result.exit_code == 0
    assert USER_EMAIL in result.stdout
