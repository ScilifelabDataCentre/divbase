"""
Tests for the "divbase-cli version" subcommand
"""

import pytest
from sqlalchemy import text
from typer.testing import CliRunner

from divbase_cli.cli_exceptions import DivBaseAPIError
from divbase_cli.divbase_cli import app

runner = CliRunner()


VERSION_1_NAME = "v1.0.0"
VERSION_2_NAME = "v2.0.0"
VERSION_3_NAME = "v3.0.0"


@pytest.fixture(autouse=True)
def clean_versions(db_session_sync):
    """
    Delete all rows from the project_versions table between each test
    """
    db_session_sync.execute(text("DELETE FROM project_version;"))
    db_session_sync.commit()

    yield


def test_add_version(logged_in_edit_user_with_existing_config):
    command = f"version add {VERSION_1_NAME}"
    result = runner.invoke(app, command)

    assert result.exit_code == 0
    assert f"New version: '{VERSION_1_NAME}'" in result.stdout


def test_add_version_with_description(logged_in_edit_user_with_existing_config):
    description = "Initial release"
    command = f'version add {VERSION_1_NAME} --description "{description}"'
    result = runner.invoke(app, command)

    assert result.exit_code == 0

    list_cmd = "version list"
    list_result = runner.invoke(app, list_cmd)

    assert VERSION_1_NAME in list_result.stdout
    assert description in list_result.stdout


def test_add_multiple_versions(logged_in_edit_user_with_existing_config):
    versions = [VERSION_1_NAME, VERSION_2_NAME, VERSION_3_NAME]
    for version in versions:
        command = f"version add {version}"
        runner.invoke(app, command)

    command = "version list"
    result = runner.invoke(app, command)

    assert result.exit_code == 0

    for version in versions:
        assert version in result.stdout


def test_attempt_add_version_that_already_exists_fails(logged_in_edit_user_with_existing_config):
    command = "version add v1.0.0"

    result = runner.invoke(app, command)
    assert result.exit_code == 0

    result = runner.invoke(app, command)
    assert result.exit_code != 0
    assert isinstance(result.exception, DivBaseAPIError)
    assert result.exception.error_type == "project_version_already_exists_error"
    assert result.exception.status_code == 400


def test_list_versions(logged_in_edit_user_with_existing_config):
    runner.invoke(app, f"version add {VERSION_1_NAME}")
    runner.invoke(app, f"version add {VERSION_2_NAME}")

    command = "version list"
    result = runner.invoke(app, command)

    assert result.exit_code == 0
    assert f"{VERSION_1_NAME}" in result.stdout
    assert f"{VERSION_2_NAME}" in result.stdout


def test_list_versions_with_and_without_include_deleted_flag(logged_in_edit_user_with_existing_config):
    """Test that deleted versions only show up when include_deleted flag is used"""
    runner.invoke(app, f"version add {VERSION_1_NAME}")
    runner.invoke(app, f"version add {VERSION_2_NAME}")
    runner.invoke(app, f"version delete {VERSION_1_NAME}")

    result = runner.invoke(app, "version list")
    assert result.exit_code == 0
    assert VERSION_1_NAME not in result.stdout
    assert VERSION_2_NAME in result.stdout

    result = runner.invoke(app, "version list --include-deleted")
    assert result.exit_code == 0
    assert VERSION_1_NAME in result.stdout
    assert VERSION_2_NAME in result.stdout


def test_list_versions_for_empty_project(logged_in_edit_user_with_existing_config):
    """Test that list version works fine if no versions exist"""
    result = runner.invoke(app, "version list")
    assert result.exit_code == 0
    assert "No versions found" in result.stdout

    result = runner.invoke(app, "version list --include-deleted")
    assert result.exit_code == 0
    assert "No versions found" in result.stdout


def test_delete_version(logged_in_edit_user_with_existing_config):
    command = f"version add {VERSION_1_NAME}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0

    command = f"version delete {VERSION_1_NAME}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert f"version: '{VERSION_1_NAME}' was deleted" in result.stdout

    list_cmd = "version list"
    result = runner.invoke(app, list_cmd)
    assert result.exit_code == 0
    assert VERSION_1_NAME not in result.stdout


def test_delete_nonexistent_version(logged_in_edit_user_with_existing_config):
    nonexistent_version = "v99.99.99"
    command = f"version delete {nonexistent_version}"
    result = runner.invoke(app, command)

    assert result.exit_code != 0
    assert isinstance(result.exception, DivBaseAPIError)
    assert result.exception.error_type == "project_version_not_found_error"
    assert result.exception.status_code == 404


def test_delete_already_deleted_version(logged_in_edit_user_with_existing_config):
    """
    Test deleting a version that's already been soft deleted
    Should not raise error, but just say it's already been deleted.
    """
    runner.invoke(app, f"version add {VERSION_1_NAME}")
    result = runner.invoke(app, f"version delete {VERSION_1_NAME}")
    assert result.exit_code == 0
    assert f"version: '{VERSION_1_NAME}' was deleted" in result.stdout

    result = runner.invoke(app, f"version delete {VERSION_1_NAME}")
    assert result.exit_code == 0
    assert VERSION_1_NAME in result.stdout
    assert "has already been soft-deleted" in result.stdout


def test_get_version_info(logged_in_edit_user_with_existing_config, CONSTANTS):
    default_project = CONSTANTS["DEFAULT_PROJECT"]
    files_in_project = CONSTANTS["PROJECT_CONTENTS"][default_project]

    command = f"version add {VERSION_1_NAME}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0

    command = f"version info {VERSION_1_NAME}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0

    assert f"{VERSION_1_NAME}" in result.stdout
    assert f"{default_project}" in result.stdout

    for file in files_in_project:
        assert f"- '{file}' :" in result.stdout


def test_get_version_info_for_version_that_does_not_exist(logged_in_edit_user_with_existing_config, CONSTANTS):
    command = "version info does_not_exist"
    result = runner.invoke(app, command)

    assert result.exit_code != 0
    assert isinstance(result.exception, DivBaseAPIError)
    assert result.exception.error_type == "project_version_not_found_error"
    assert result.exception.status_code == 404


def test_get_version_info_for_deleted_version(logged_in_edit_user_with_existing_config, CONSTANTS):
    """Should still work"""
    default_project = CONSTANTS["DEFAULT_PROJECT"]
    files_in_project = CONSTANTS["PROJECT_CONTENTS"][default_project]

    runner.invoke(app, f"version add {VERSION_1_NAME}")
    runner.invoke(app, f"version delete {VERSION_1_NAME}")

    command = f"version info {VERSION_1_NAME}"
    result = runner.invoke(app, command)

    assert result.exit_code == 0
    assert VERSION_1_NAME in result.stdout
    assert "WARNING: This version has been soft-deleted" in result.stdout

    for file in files_in_project:
        assert file in result.stdout


def test_get_version_updates_hashes_on_new_upload(logged_in_edit_user_with_existing_config, CONSTANTS, fixtures_dir):
    """
    Test a simple protocol where:
    1. upload file
    2. add a version
    3. upload the same file again (disabling the safe mode setting so it can be re-uploaded)
    4. add another version

    Validate that the hashes for the uploaded file in each version are different.
    """
    test_file_name = CONSTANTS["FILES_TO_UPLOAD_DOWNLOAD"][0]
    test_file_path = (fixtures_dir / test_file_name).resolve()
    upload_cmd = f"files upload {test_file_path} --disable-safe-mode"

    result = runner.invoke(app, f"version add {VERSION_3_NAME}")
    assert result.exit_code == 0

    result = runner.invoke(app, upload_cmd)
    assert result.exit_code == 0

    result = runner.invoke(app, f"version add {VERSION_1_NAME}")
    assert result.exit_code == 0

    result = runner.invoke(app, upload_cmd)
    assert result.exit_code == 0

    result = runner.invoke(app, f"version add {VERSION_2_NAME}")
    assert result.exit_code == 0

    info_result_1 = runner.invoke(app, f"version info {VERSION_1_NAME}")
    assert info_result_1.exit_code == 0
    info_result_2 = runner.invoke(app, f"version info {VERSION_2_NAME}")
    assert info_result_2.exit_code == 0

    hash_v1, hash_v2 = "", ""
    for line in info_result_1.stdout.splitlines():
        if test_file_name in line:
            hash_v1 = line.split(f"- '{test_file_name}' : ")[1].strip()
            break
    for line in info_result_2.stdout.splitlines():
        if test_file_name in line:
            hash_v2 = line.split(f"- '{test_file_name}' : ")[1].strip()
            break

    assert hash_v1 != "", "Hash for version 1 not found"
    assert hash_v2 != "", "Hash for version 2 not found"
    assert hash_v1 != hash_v2, "Hashes for the same file in different versions should be different"
