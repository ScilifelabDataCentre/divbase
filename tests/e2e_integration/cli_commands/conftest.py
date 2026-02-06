"""
Setup for pytest fixtures for e2e testing of the CLI commands.

Pytest fixtures set up (and tear down) a test environment with the full DivBase stack running locally,

The S3FileManager class is patched in all tests to use this test MinIO server,
it is autoused, so it does not need to be specified in each test.
"""

import logging
from pathlib import Path

import pytest
from typer.testing import CliRunner

from divbase_cli.cli_config import cli_settings
from divbase_cli.divbase_cli import app

runner = CliRunner()

logger = logging.getLogger(__name__)


@pytest.fixture
def tmp_config_path(tmp_path):
    """
    Fixture to provide a path to where the a configuration can be created.
    """
    return tmp_path / "test_config.yaml"


@pytest.fixture
def logged_out_user_with_fresh_config():
    """
    Fixture to provide a user with a fresh config file that is not logged in anywhere.
    """
    # ensure no config or tokens file exist before test
    cli_settings.CONFIG_PATH.unlink(missing_ok=True)
    cli_settings.TOKENS_PATH.unlink(missing_ok=True)
    yield
    # clean up after test, delete config and tokens file
    cli_settings.CONFIG_PATH.unlink(missing_ok=True)
    cli_settings.TOKENS_PATH.unlink(missing_ok=True)


@pytest.fixture
def logged_out_user_with_existing_config(CONSTANTS):
    """
    Fixture to provide a NOT logged in user with an "existing" user configuration file with
    some existing projects and a default project set.
    """
    # ensure no config or tokens file exist before test
    cli_settings.CONFIG_PATH.unlink(missing_ok=True)
    cli_settings.TOKENS_PATH.unlink(missing_ok=True)

    # running any cmd that requires the config file will create it
    # so we can just run the config add cmd directly.
    for project in CONSTANTS["PROJECT_TO_BUCKET_MAP"]:
        add_command = f"config add {project}"
        result = runner.invoke(app, add_command)
        assert result.exit_code == 0

    set_default_command = f"config set-default {CONSTANTS['DEFAULT_PROJECT']}"
    result = runner.invoke(app, set_default_command)
    assert result.exit_code == 0

    yield

    # clean up after test, delete config and tokens file
    cli_settings.CONFIG_PATH.unlink(missing_ok=True)
    cli_settings.TOKENS_PATH.unlink(missing_ok=True)


@pytest.fixture
def logged_in_admin_with_existing_config(CONSTANTS):
    """Fixture to provide a logged in admin user with existing config."""
    yield from _create_logged_in_user_fixture("admin")(CONSTANTS)


@pytest.fixture
def logged_in_read_user_with_existing_config(CONSTANTS):
    """Fixture to provide a logged in read user with existing config."""
    yield from _create_logged_in_user_fixture("read user")(CONSTANTS)


@pytest.fixture
def logged_in_edit_user_with_existing_config(CONSTANTS):
    """Fixture to provide a logged in edit user with existing config."""
    yield from _create_logged_in_user_fixture("edit user")(CONSTANTS)


@pytest.fixture
def logged_in_manage_user_with_existing_config(CONSTANTS):
    """Fixture to provide a logged in manage user with existing config."""
    yield from _create_logged_in_user_fixture("manage user")(CONSTANTS)


@pytest.fixture
def logged_in_edit_user_query_project_only_with_existing_config(CONSTANTS):
    """Fixture to provide a logged in edit user (who only belongs to query-project) with existing config."""
    yield from _create_logged_in_user_fixture("edit user query-project only")(CONSTANTS)


@pytest.fixture
def logged_in_manage_user_query_project_only_with_existing_config(CONSTANTS):
    """Fixture to provide a logged in manage user (who only belongs to query-project) with existing config."""
    yield from _create_logged_in_user_fixture("manage user query-project only")(CONSTANTS)


def _create_logged_in_user_fixture(user_type: str):
    """
    Factory function to create a logged-in user fixture for a specific user type.

    Args:
        user_type: One of "admin", "read user", "edit user", "manage user"
    """

    def factory(CONSTANTS):
        # ensure no config or tokens file exist before test
        cli_settings.CONFIG_PATH.unlink(missing_ok=True)
        cli_settings.TOKENS_PATH.unlink(missing_ok=True)

        # running any cmd that requires the config file will create it
        for project in CONSTANTS["PROJECT_TO_BUCKET_MAP"]:
            add_command = f"config add {project}"
            result = runner.invoke(app, add_command)
            assert result.exit_code == 0

        set_default_command = f"config set-default {CONSTANTS['DEFAULT_PROJECT']}"
        result = runner.invoke(app, set_default_command)
        assert result.exit_code == 0

        # Get credentials based on user type
        if user_type == "admin":
            user_creds = CONSTANTS["ADMIN_CREDENTIALS"]
        else:
            user_creds = CONSTANTS["TEST_USERS"][user_type]

        login_command = f"auth login {user_creds['email']} --password {user_creds['password']}"
        result = runner.invoke(app, login_command)
        assert result.exit_code == 0, f"Login failed: {result.output}"

        yield

        # clean up after test, delete config and tokens file
        cli_settings.CONFIG_PATH.unlink(missing_ok=True)
        cli_settings.TOKENS_PATH.unlink(missing_ok=True)

    return factory


@pytest.fixture
def fixtures_dir():
    """Path to the fixtures directory."""
    return Path(__file__).parent.parent.parent / "fixtures"
