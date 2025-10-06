"""
E2E tests for the "divbase-cli auth" CLI commands.

NOTE: Without a user config file already created, these tests will fail

TODO - add some non e2e tests for logging in with expired tokens etc.
"""

import pytest
from typer.testing import CliRunner

from divbase_cli.divbase_cli import app
from divbase_lib.exceptions import AuthenticationError, DivBaseAPIConnectionError

runner = CliRunner()

# TODO
# move these over to CONSTANTS later..
# Create non admin user for testing
admin_credentials = {
    "email": "admin@divbase.com",
    "password": "badpassword",
}


@pytest.fixture(autouse=True)
def ensure_logged_out():
    """Fixture to ensure we are logged out before and after each test."""
    command = "auth logout"
    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert "Logged out successfully" in result.stdout
    yield
    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert "Logged out successfully" in result.stdout


def log_in_as_admin():
    """
    Helper function to log in as admin user.
    Not a fixture as want to ensure logged out before and after each test and could be timing issues when
    combined with other fixtures (e.g. ensure_logged_out).
    """
    command = f"auth login {admin_credentials['email']} --password {admin_credentials['password']}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert "Logged in successfully" in result.stdout


def test_login_command(fresh_config_path):
    command = f"auth login {admin_credentials['email']} --password {admin_credentials['password']}"

    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert "Logged in successfully" in result.stdout
    assert admin_credentials["email"] in result.stdout


def test_login_command_with_password_prompted(fresh_config_path):
    command = f"auth login {admin_credentials['email']}"

    result = runner.invoke(app, command, input=f"{admin_credentials['password']}\n")
    assert result.exit_code == 0
    assert "Logged in successfully" in result.stdout
    assert admin_credentials["email"] in result.stdout


def test_login_command_fails_with_invalid_credentials(fresh_config_path):
    """Test login command fails with invalid credentials."""
    command = f"auth login {admin_credentials['email']} --password wrongpassword"

    result = runner.invoke(app, command)
    assert result.exit_code != 0
    assert isinstance(result.exception, AuthenticationError)
    assert "Invalid email or password" in str(result.exception)


def test_login_command_with_invalid_server_url(fresh_config_path):
    """Test login command fails with an invalid server URL."""
    command = f"auth login {admin_credentials['email']} --password {admin_credentials['password']} --divbase-url https://invalid-url"

    result = runner.invoke(app, command)
    assert result.exit_code != 0
    assert isinstance(result.exception, DivBaseAPIConnectionError)


def test_login_command_already_logged_in(fresh_config_path):
    """Test login command when already logged in."""
    log_in_as_admin()

    command = f"auth login {admin_credentials['email']} --password {admin_credentials['password']}"

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
    assert admin_credentials["email"] in result.stdout


def test_force_login_option(fresh_config_path):
    """Should not prompt about logging in again"""
    log_in_as_admin()
    command = f"auth login {admin_credentials['email']} --password {admin_credentials['password']} --force"

    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert "Logged in successfully" in result.stdout
    assert admin_credentials["email"] in result.stdout


def test_logout_command(fresh_config_path):
    """Test basic usage of logout and that running multiple times is ok."""
    command = "auth logout"

    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert "Logged out successfully" in result.stdout

    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert "Logged out successfully" in result.stdout


def test_login_logout_cycle(fresh_config_path):
    """Test a few repeated login/logout cycles."""
    login_command = f"auth login {admin_credentials['email']} --password {admin_credentials['password']}"
    logout_command = "auth logout"

    for _ in range(3):
        result = runner.invoke(app, login_command)
        assert result.exit_code == 0
        assert "Logged in successfully" in result.stdout
        assert admin_credentials["email"] in result.stdout

        result = runner.invoke(app, logout_command)
        assert result.exit_code == 0
        assert "Logged out successfully" in result.stdout


def test_whoami_command(fresh_config_path):
    """Test basic usage of whoami command."""
    log_in_as_admin()
    command = "auth whoami"

    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert admin_credentials["email"] in result.stdout


def test_whoami_command_fails_if_not_logged_in(fresh_config_path):
    """Test basic usage of whoami command."""
    command = "auth whoami"

    result = runner.invoke(app, command)
    assert result.exit_code != 0
    assert isinstance(result.exception, AuthenticationError)
