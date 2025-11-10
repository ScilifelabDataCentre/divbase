"""
Tests for the "divbase-cli files" commands

All tests are run against the docker compose test overlay

A clean project (its bucket is auto emptied before each test) is available to any test that requires a clean slate.
"""

import boto3
import pytest
from typer.testing import CliRunner

from divbase_cli.divbase_cli import app
from divbase_lib.exceptions import FilesAlreadyInBucketError

runner = CliRunner()


@pytest.fixture(autouse=True)
def start_with_clean_project(CONSTANTS):
    """
    For tests that require a project with a clean bucket, this fixture will
    ensure that the CONSTANTS["CLEANED_PROJECT"]'s bucket is empty before running the test.

    Caution:
    If you modify the approach make sure your implementation does not just add delete markers.
    The files need to be actually deleted.
    """
    s3_resource = boto3.resource(
        "s3",
        endpoint_url=CONSTANTS["MINIO_URL"],
        aws_access_key_id=CONSTANTS["BAD_ACCESS_KEY"],
        aws_secret_access_key=CONSTANTS["BAD_SECRET_KEY"],
    )

    # pylance does not understand boto3 resource returns types, hence ignore below
    bucket = s3_resource.Bucket(CONSTANTS["CLEANED_PROJECT"])  # type: ignore
    bucket.object_versions.delete()

    yield


def test_list_files(logged_in_edit_user_with_existing_config, CONSTANTS):
    """Test basic usage of files list command."""
    command = "files list"

    result = runner.invoke(app, command)
    assert result.exit_code == 0

    default_project = CONSTANTS["DEFAULT_PROJECT"]

    for file in CONSTANTS["PROJECT_CONTENTS"][default_project]:
        assert file in result.stdout, f"File {file} not found in the output of the list_files command"


def test_list_non_default_project(logged_in_edit_user_with_existing_config, CONSTANTS):
    """Test list files for the non-default project."""
    non_default_project = CONSTANTS["NON_DEFAULT_PROJECT"]
    files_in_project = CONSTANTS["PROJECT_CONTENTS"][non_default_project]

    command = f"files list --project {non_default_project}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0
    for file in files_in_project:
        assert file in result.stdout


def test_list_files_empty_project(logged_in_edit_user_with_existing_config, CONSTANTS):
    """Test list files for an empty project."""
    command = f"files list --project {CONSTANTS['EMPTY_PROJECT']}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert "No files found" in result.stdout


def test_upload_1_file(logged_in_edit_user_with_existing_config, tmp_path):
    """Test upload 1 file to the project."""
    test_file = tmp_path / "fake_test_file.txt"
    test_file.write_text("testing, testing 1 2 3...")

    command = f"files upload {test_file}"
    result = runner.invoke(app, command)

    assert result.exit_code == 0
    assert f"{str(test_file)}" in result.stdout


def test_upload_1_file_to_non_default_project(logged_in_edit_user_with_existing_config, CONSTANTS, fixtures_dir):
    """Specify a project when uploading a file."""
    test_file = (fixtures_dir / CONSTANTS["FILES_TO_UPLOAD_DOWNLOAD"][0]).resolve()

    command = f"files upload {test_file} --project {CONSTANTS['CLEANED_PROJECT']}"
    result = runner.invoke(app, command)

    assert result.exit_code == 0
    assert f"{str(test_file)}" in result.stdout


def test_upload_multiple_files_at_once(logged_in_edit_user_with_existing_config, CONSTANTS, fixtures_dir):
    test_files = [(fixtures_dir / file_name).resolve() for file_name in CONSTANTS["FILES_TO_UPLOAD_DOWNLOAD"]]

    command = f"files upload {' '.join(map(str, test_files))} --project {CONSTANTS['CLEANED_PROJECT']}"
    result = runner.invoke(app, command)

    assert result.exit_code == 0
    for file in test_files:
        assert f"{str(file)}" in result.stdout


def test_upload_dir_contents(logged_in_edit_user_with_existing_config, CONSTANTS, fixtures_dir):
    """Test upload all files in a directory."""
    files = [x for x in fixtures_dir.glob("*") if x.is_file()]  # does not get subdirs

    command = f"files upload --upload-dir {fixtures_dir.resolve()} --project {CONSTANTS['CLEANED_PROJECT']}"
    result = runner.invoke(app, command)

    assert result.exit_code == 0
    clean_stdout = result.stdout.replace("\n", "")  # newlines can cause issues in the assert below
    for file in files:
        assert str(file.resolve()) in clean_stdout


def test_reupload_same_file_fails(logged_in_edit_user_with_existing_config, CONSTANTS, fixtures_dir):
    """Test upload with safe mode on (default) works the first time, but fails on subsequent attempts."""
    file_name = CONSTANTS["FILES_TO_UPLOAD_DOWNLOAD"][0]
    file_path = f"{fixtures_dir}/{file_name}"
    command = f"files upload {file_path} --project {CONSTANTS['CLEANED_PROJECT']}"

    result = runner.invoke(app, command)
    assert result.exit_code == 0

    result = runner.invoke(app, command)
    assert result.exit_code != 0
    assert isinstance(result.exception, FilesAlreadyInBucketError)


def test_reupload_of_same_file_with_safe_mode_off_works(
    logged_in_edit_user_with_existing_config, CONSTANTS, fixtures_dir
):
    """Test upload with safe mode off works for reuploading the same file."""
    file_name = CONSTANTS["FILES_TO_UPLOAD_DOWNLOAD"][0]
    file_path = f"{fixtures_dir}/{file_name}"
    command = f"files upload {file_path} --project {CONSTANTS['CLEANED_PROJECT']} --disable-safe-mode"

    result = runner.invoke(app, command)
    assert result.exit_code == 0

    result = runner.invoke(app, command)
    assert result.exit_code == 0


def test_no_file_uploaded_if_some_duplicated(logged_in_edit_user_with_existing_config, CONSTANTS, fixtures_dir):
    """
    Test that no files are uploaded with safe mode on (which is the default)
    if at least 1 of the files trying to be uploaded already exists in the project's bucket.
    """
    test_files = [(fixtures_dir / file_name).resolve() for file_name in CONSTANTS["FILES_TO_UPLOAD_DOWNLOAD"]]

    # upload just 1 of the files first
    command = f"files upload {test_files[0]} --project {CONSTANTS['CLEANED_PROJECT']}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0

    # none should be uploaded as the first one already exists
    command = f"files upload {' '.join(map(str, test_files))} --project {CONSTANTS['CLEANED_PROJECT']}"
    result = runner.invoke(app, command)
    assert result.exit_code != 0
    assert isinstance(result.exception, FilesAlreadyInBucketError)

    command = f"files list --project {CONSTANTS['CLEANED_PROJECT']}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert test_files[0].name in result.stdout
    for file in test_files[1:]:
        assert file.name not in result.stdout, f"File {file.name} was uploaded when it shouldn't have been."


def test_download_1_file(logged_in_edit_user_with_existing_config, CONSTANTS, tmp_path):
    file_name = CONSTANTS["PROJECT_CONTENTS"][CONSTANTS["DEFAULT_PROJECT"]][0]
    download_dir = tmp_path / "downloads"
    download_dir.mkdir()

    command = f"files download {file_name} --download-dir {download_dir}"
    result = runner.invoke(app, command)

    assert result.exit_code == 0
    assert file_name in result.stdout
    assert (download_dir / file_name).exists()


def test_download_multiple_files(logged_in_edit_user_with_existing_config, CONSTANTS, tmp_path):
    files_in_project = CONSTANTS["PROJECT_CONTENTS"][CONSTANTS["DEFAULT_PROJECT"]]
    download_dir = tmp_path / "downloads"
    download_dir.mkdir()

    command = f"files download {' '.join(files_in_project)} --download-dir {download_dir}"
    result = runner.invoke(app, command)

    assert result.exit_code == 0
    for file_name in files_in_project:
        assert file_name in result.stdout
        assert (download_dir / file_name).exists()


def test_download_from_non_default_project(logged_in_edit_user_with_existing_config, CONSTANTS, tmp_path):
    non_default_project = CONSTANTS["NON_DEFAULT_PROJECT"]
    file_to_download = CONSTANTS["PROJECT_CONTENTS"][non_default_project][0]

    download_dir = tmp_path / "downloads"
    download_dir.mkdir()

    command = f"files download {file_to_download} --project {non_default_project} --download-dir {download_dir}"
    result = runner.invoke(app, command)

    assert result.exit_code == 0
    assert file_to_download in result.stdout
    assert (download_dir / file_to_download).exists()


def test_download_using_file_list(logged_in_edit_user_with_existing_config, CONSTANTS, tmp_path):
    files_in_project = CONSTANTS["PROJECT_CONTENTS"][CONSTANTS["DEFAULT_PROJECT"]]
    download_dir = tmp_path / "downloads"
    download_dir.mkdir()

    file_list = tmp_path / "file_list.txt"
    with open(file_list, "w") as f:
        f.write("\n".join(files_in_project))

    command = f"files download --file-list {file_list} --download-dir {download_dir}"
    result = runner.invoke(app, command)

    assert result.exit_code == 0
    for file_name in files_in_project:
        assert file_name in result.stdout
        assert (download_dir / file_name).exists()


def test_download_nonexistent_file(logged_in_edit_user_with_existing_config, tmp_path):
    download_dir = tmp_path / "downloads"
    download_dir.mkdir()

    command = f"files download nonexistent_file.txt --download-dir {download_dir}"
    result = runner.invoke(app, command)

    assert result.exit_code == 0
    assert "WARNING: Failed to download the following files:" in result.stdout
    assert "nonexistent_file.txt" in result.stdout


def test_download_at_a_project_version(logged_in_edit_user_with_existing_config, CONSTANTS, tmp_path, fixtures_dir):
    """Test downloading at specified project versions."""
    clean_project = CONSTANTS["CLEANED_PROJECT"]

    download_dir = tmp_path / "download_v1"
    download_dir.mkdir()

    file_name = "test_file.txt"
    file_path = tmp_path / file_name
    v1_content = "This is version 1.0.0 content"
    v2_content = "This is version 2.0.0 content - UPDATED"

    # Create file, upload v1 contents and create project version v1.0.0
    with open(file_path, "w") as f:
        f.write(v1_content)

    command = f"files upload {file_path} --project {clean_project}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0

    command = f"version add v1.0.0 --project {clean_project}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0

    # Create file, upload v2 contents and create project version v2.0.0
    with open(file_path, "w") as f:
        f.write(v2_content)

    command = f"files upload {file_path} --project {clean_project}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0

    command = f"version add v2.0.0 --project {clean_project}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0

    # Download + assert contents at v1.0.0
    command = (
        f"files download {file_name} --bucket-version v1.0.0 --project {clean_project} --download-dir {download_dir}"
    )
    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert (download_dir / file_name).exists()

    with open(download_dir / file_name) as f:
        downloaded_content = f.read()
    assert downloaded_content == v1_content

    # Download + assert contents at v2.0.0
    command = (
        f"files download {file_name} --bucket-version v2.0.0 --project {clean_project} --download-dir {download_dir}"
    )
    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert (download_dir / file_name).exists()

    with open(download_dir / file_name) as f:
        downloaded_content = f.read()
    assert downloaded_content == v2_content

    # Download without specifiying bucket version works like "latest"
    command = f"files download {file_name} --project {clean_project} --download-dir {download_dir}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert (download_dir / file_name).exists()

    with open(download_dir / file_name) as f:
        downloaded_content = f.read()
    assert downloaded_content == v2_content


def test_remove_with_dry_run(logged_in_edit_user_with_existing_config, CONSTANTS):
    file_name = CONSTANTS["PROJECT_CONTENTS"][CONSTANTS["DEFAULT_PROJECT"]][0]

    command = f"files remove {file_name} --dry-run"
    result = runner.invoke(app, command)

    assert result.exit_code == 0
    assert f"{file_name}" in result.stdout

    command = "files list"
    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert file_name in result.stdout


def test_remove_file(logged_in_edit_user_with_existing_config, CONSTANTS, fixtures_dir):
    """Test removing a file from the project's bucket, using a clean project to avoid side effects in other tests."""
    clean_project = CONSTANTS["CLEANED_PROJECT"]

    file_name = CONSTANTS["FILES_TO_UPLOAD_DOWNLOAD"][0]
    file_path = f"{fixtures_dir}/{file_name}"

    command = f"files upload {file_path} --project {clean_project}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0

    command = f"files list --project {clean_project}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert file_name in result.stdout

    command = f"files remove {file_name} --project {clean_project}"
    result = runner.invoke(app, command)

    assert result.exit_code == 0
    assert f"{file_name}" in result.stdout

    command = f"files list --project {clean_project}"
    result = runner.invoke(app, command)
    assert result.exit_code == 0
    assert file_name not in result.stdout
