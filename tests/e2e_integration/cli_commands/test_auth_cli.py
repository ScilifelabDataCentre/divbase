"""
E2E tests for the "divbase-cli auth" CLI commands.

NOTE: Without a user config file already created, these tests will fail
"""

import shutil

from typer.testing import CliRunner

from divbase_cli.cli_config import cli_settings
from divbase_cli.cli_exceptions import AuthenticationError, DivBaseAPIConnectionError
from divbase_cli.divbase_cli import app
from divbase_cli.user_auth import LOGIN_AGAIN_MESSAGE

runner = CliRunner()

USER_EMAIL = "edit@divbase.se"
USER_PASSWORD = "badpassword"


def log_in_as_user():
    """
    Helper function to log in as a user.
    Not a fixture as want to ensure logged out before and after each test and could be timing issues when
    combined with other fixtures (e.g. ensure_logged_out).
    """
    command = f"auth login {USER_EMAIL} --password {USER_PASSWORD}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert "Logged in successfully" in result.stdout


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


def test_login_command(logged_out_user_with_fresh_config):
    command = f"auth login {USER_EMAIL} --password {USER_PASSWORD}"

    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert "Logged in successfully" in result.stdout
    assert USER_EMAIL in result.stdout


def test_login_command_with_password_prompted(logged_out_user_with_fresh_config):
    command = f"auth login {USER_EMAIL}"

    result = runner.invoke(app, command, input=f"{USER_PASSWORD}\n")
    assert result.exit_code == 0
    assert "Logged in successfully" in result.stdout
    assert USER_EMAIL in result.stdout


def test_login_command_fails_with_invalid_credentials(logged_out_user_with_fresh_config):
    """Test login command fails with invalid credentials."""
    command = f"auth login {USER_EMAIL} --password wrongpassword"

    result = runner.invoke(app, command)
    assert result.exit_code != 0
    assert isinstance(result.exception, AuthenticationError)
    assert "Invalid email or password" in str(result.exception)


def test_login_command_with_invalid_server_url(logged_out_user_with_fresh_config):
    """Test login command fails with an invalid server URL."""
    command = f"auth login {USER_EMAIL} --password {USER_PASSWORD} --divbase-url https://invalid-url"

    result = runner.invoke(app, command)
    assert result.exit_code != 0
    assert isinstance(result.exception, DivBaseAPIConnectionError)


def test_login_command_already_logged_in(logged_out_user_with_fresh_config):
    """Test login command when already logged in."""
    log_in_as_user()

    command = f"auth login {USER_EMAIL} --password {USER_PASSWORD}"

    # cancel login, when warned already logged in
    result = runner.invoke(app, command, input="n\n")
    assert result.exit_code == 0
    assert "Already logged in to" in result.stdout
    assert "Do you want to login again?" in result.stdout
    assert "Login cancelled." in result.stdout

    # don't cancel login process when warned already logged in
    result = runner.invoke(app, command, input="y\n")
    assert result.exit_code == 0
    assert "Already logged in to" in result.stdout
    assert "Do you want to login again?" in result.stdout
    assert "Logged in successfully" in result.stdout
    assert USER_EMAIL in result.stdout


def test_force_login_option(logged_out_user_with_fresh_config):
    """Should not prompt about logging in again"""
    log_in_as_user()
    command = f"auth login {USER_EMAIL} --password {USER_PASSWORD} --force"

    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert "Logged in successfully" in result.stdout
    assert USER_EMAIL in result.stdout


def test_logout_command(logged_out_user_with_fresh_config):
    """Test basic usage of logout and that running multiple times is ok."""
    command = "auth logout"

    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert "Logged out successfully" in result.stdout

    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert "Logged out successfully" in result.stdout


def test_login_logout_cycle(logged_out_user_with_fresh_config):
    """Test a few repeated login/logout cycles."""
    login_command = f"auth login {USER_EMAIL} --password {USER_PASSWORD}"
    logout_command = "auth logout"

    for _ in range(3):
        result = runner.invoke(app, login_command)
        assert result.exit_code == 0
        assert "Logged in successfully" in result.stdout
        assert USER_EMAIL in result.stdout

        result = runner.invoke(app, logout_command)
        assert result.exit_code == 0
        assert "Logged out successfully" in result.stdout


def test_whoami_command(logged_out_user_with_fresh_config):
    """Test basic usage of whoami command."""
    log_in_as_user()
    command = "auth whoami"

    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert USER_EMAIL in result.stdout


def test_whoami_command_fails_if_not_logged_in(logged_out_user_with_fresh_config):
    """Test basic usage of whoami command."""
    command = "auth whoami"

    result = runner.invoke(app, command)
    assert result.exit_code != 0
    assert isinstance(result.exception, AuthenticationError)


def test_whoami_command_needing_refresh_token(logged_out_user_with_fresh_config):
    """
    We simulate that the access token has expired by manually setting it to be expired in the users .secrets file

    This will force the CLI to use the refresh token to get a new access token when running the whoami command.
    """
    log_in_as_user()
    make_tokens_expired(access=True)

    whoami_command = "auth whoami"
    result = runner.invoke(app, whoami_command)
    assert result.exit_code == 0
    assert USER_EMAIL in result.stdout


def test_whoami_with_expired_tokens_fails(logged_in_admin_with_existing_config):
    """
    Simulate that both the access token and refresh token have expired by manually setting them to be expired
    in the users .secrets file.

    This should force authentication to fail when running the whoami command, requiring the user to log in again.
    """
    command = "auth whoami"

    make_tokens_expired(access=True, refresh=True)

    result = runner.invoke(app, command)
    assert result.exit_code != 0
    assert isinstance(result.exception, AuthenticationError)
    assert LOGIN_AGAIN_MESSAGE in str(result.exception)


def test_using_revoked_refresh_token_fails(logged_out_user_with_fresh_config, tmp_path):
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

    result = runner.invoke(app, logout_command)
    assert result.exit_code == 0

    shutil.copy(tmp_path / "tokens_backup.yaml", cli_settings.TOKENS_PATH)
    shutil.copy(tmp_path / "config_backup.yaml", cli_settings.CONFIG_PATH)
    make_tokens_expired(access=True)

    result = runner.invoke(app, whoami_command)
    assert result.exit_code != 0
    assert isinstance(result.exception, AuthenticationError)
    assert LOGIN_AGAIN_MESSAGE in str(result.exception)

    log_in_as_user()

    result = runner.invoke(app, whoami_command)
    assert result.exit_code == 0
    assert USER_EMAIL in result.stdout
