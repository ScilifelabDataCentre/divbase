"""
Setup for pytest fixtures for e2e testing of the CLI commands.

Pytest fixtures set up (and tear down) a test environment with the full DivBase stack running locally,

The S3FileManager class is patched in all tests to use this test MinIO server,
it is autoused, so it does not need to be specified in each test.
"""

import contextlib
import logging

import boto3
import keyring
import pytest
from keyring.errors import KeyringError
from typer.testing import CliRunner

from divbase_cli.cli_config import cli_settings
from divbase_cli.divbase_cli import app

runner = CliRunner()

logger = logging.getLogger(__name__)


@pytest.fixture(autouse=True)
def no_pat(monkeypatch):
    """
    Ensures DIVBASE_API_PAT is never taken from test runner environment (aka dev's token),
    Could interfere otherwise.

    Specific tests can instead override this fixture if they want to test behavior with a PAT.
    """
    monkeypatch.setattr(cli_settings, "DIVBASE_API_PAT", None)


@pytest.fixture
def tmp_config_path(tmp_path):
    """
    Fixture to provide a path to where the a configuration can be created.
    """
    return tmp_path / "test_config.yaml"


@pytest.fixture
def logged_out_user_with_existing_config(CONSTANTS):
    """
    Fixture to provide a NOT logged in user with an "existing" user configuration file with
    some existing projects and a default project set.

    NOTE: config and tokens cleaned up before and after tests by clean_tmp_config_and_tokens_between_tests fixture
    """
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


@pytest.fixture
def logged_in_admin_with_existing_config(CONSTANTS):
    """Fixture to provide a logged in admin user with existing config."""
    yield from _create_logged_in_user_fixture("admin")(CONSTANTS)


@pytest.fixture
def logged_in_read_user_with_existing_config(CONSTANTS):
    """Fixture to provide a logged in read user with existing config."""
    yield from _create_logged_in_user_fixture("read user")(CONSTANTS)


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

    NOTE: Whilst config and tokens cleaned up before and after tests by clean_tmp_config_and_tokens_between_tests fixture
    This factory can be used multiple times in a single test, hence need to clean up existing config and tokens here too.
    """

    def factory(CONSTANTS):
        # ensure no config or tokens file exist before test
        cli_settings.CONFIG_PATH.unlink(missing_ok=True)
        # tokens can either be stored in device keyring (or in a fallback file if e.g. keyring not available - likely for CI or disabled for a test)
        with contextlib.suppress(KeyringError):
            keyring.delete_password(service_name=cli_settings.KEYRING_SERVICE, username=cli_settings.KEYRING_USERNAME)
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

        login_command = f"auth login {user_creds['email']}"
        result = runner.invoke(app=app, args=login_command, input=f"{user_creds['password']}\n")
        assert result.exit_code == 0, f"Login failed: {result.output}"

        yield

        # clean up after test, delete config and tokens file
        cli_settings.CONFIG_PATH.unlink(missing_ok=True)
        # tokens can either be stored in device keyring (or in a fallback file if e.g. keyring not available - likely for CI or disabled for a test)
        with contextlib.suppress(KeyringError):
            keyring.delete_password(service_name=cli_settings.KEYRING_SERVICE, username=cli_settings.KEYRING_USERNAME)
        cli_settings.TOKENS_PATH.unlink(missing_ok=True)

    return factory


@pytest.fixture
def cleaned_project_bucket(CONSTANTS):
    """
    Ensure cleaned-project bucket is empty before and after test.

    Use this for tests that require deterministic project bucket contents.
    """
    s3_resource = boto3.resource(
        "s3",
        endpoint_url=CONSTANTS["MINIO_URL"],
        aws_access_key_id=CONSTANTS["BAD_ACCESS_KEY"],
        aws_secret_access_key=CONSTANTS["BAD_SECRET_KEY"],
    )

    project_name = CONSTANTS["CLEANED_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    # pylance does not understand boto3 resource return types.
    bucket = s3_resource.Bucket(bucket_name)  # type: ignore
    bucket.object_versions.delete()

    yield project_name, bucket_name

    bucket.object_versions.delete()
