"""Tests for the cleanup_expired_files standalone script."""

import contextlib
from datetime import datetime, timedelta, timezone
from pathlib import Path
from unittest.mock import patch

import pytest
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session

from divbase_api.models.project_versions import ProjectVersionDB
from divbase_api.scripts import cleanup_expired_files
from divbase_api.scripts.cleanup_expired_files import (
    DEFAULT_RESULTS_FILES_RETENTION_DAYS,
    DEFAULT_SOFT_DELETED_FILES_RETENTION_DAYS,
)
from divbase_api.services.s3_client import S3FileManager, create_s3_file_manager
from divbase_lib.divbase_constants import QUERY_RESULTS_FILE_PREFIX


@pytest.fixture(autouse=True)
def start_with_clean_project(cleaned_project_bucket):
    """
    For tests that require a project with a clean bucket, this fixture will
    ensure that the CONSTANTS["CLEANED_PROJECT"]'s bucket is empty before and after running the test.
    """
    yield


@pytest.fixture(scope="module")
def s3_file_manager(CONSTANTS) -> S3FileManager:
    """Provides an S3FileManager instance configured for the test environment."""
    return create_s3_file_manager(url=CONSTANTS["MINIO_URL"])


def _run_script_with_mocked_time(mocked_now: datetime) -> None:
    """Run cleanup_expired_files.main() with datetime mocked to mocked_now."""
    with (
        patch("divbase_api.scripts.cleanup_expired_files.datetime") as mock_dt_script,
        patch("divbase_api.services.s3_client.datetime") as mock_dt_s3,
    ):
        mock_dt_script.now.return_value = mocked_now
        mock_dt_s3.now.return_value = mocked_now
        cleanup_expired_files.main()


@pytest.mark.parametrize(
    "days_offset,expect_deleted",
    [
        (DEFAULT_SOFT_DELETED_FILES_RETENTION_DAYS + 1, True),
        (DEFAULT_SOFT_DELETED_FILES_RETENTION_DAYS - 1, False),
    ],
)
def test_hard_delete_expired_files_script(
    db_session_sync: Session,
    s3_file_manager: S3FileManager,
    project_map: dict,
    CONSTANTS: dict,
    tmp_path: Path,
    days_offset: int,
    expect_deleted: bool,
):
    """
    Test that the hard_delete_expired_files script:
        - Hard-deletes expired soft-deleted objects (deletion marker older than cutoff).
        - Keeps S3 versions that are protected by a ProjectVersionDB entry.
        - Re-adds a deletion marker when a protected version would otherwise become the latest.
    """
    project_name = CONSTANTS["CLEANED_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    project_id = project_map[project_name]

    # Upload two versions of a file that will be fully purged (no protected versions).
    file_to_purge = tmp_path / "test_purge_full.txt"
    file_to_purge.write_text("delete me v1")
    s3_file_manager.upload_files(bucket_name=bucket_name, to_upload={file_to_purge.name: file_to_purge})
    file_to_purge.write_text("delete me v2")
    s3_file_manager.upload_files(bucket_name=bucket_name, to_upload={file_to_purge.name: file_to_purge})

    # Upload two versions of a file where v1 will be protected by a project version.
    file_to_keep = tmp_path / "test_purge_protected.txt"
    file_to_keep.write_text("in a project version, so kept")
    s3_file_manager.upload_files(bucket_name=bucket_name, to_upload={file_to_keep.name: file_to_keep})
    v1_details = s3_file_manager.state_of_latest_version_of_all_files(bucket_name=bucket_name)[file_to_keep.name]

    file_to_keep.write_text("will be deleted but prior version will be kept")
    s3_file_manager.upload_files(bucket_name=bucket_name, to_upload={file_to_keep.name: file_to_keep})
    v2_details = s3_file_manager.state_of_latest_version_of_all_files(bucket_name=bucket_name)[file_to_keep.name]

    s3_file_manager.soft_delete_objects(objects=[file_to_keep.name, file_to_purge.name], bucket_name=bucket_name)

    # Protect v1 of file_to_keep via a ProjectVersionDB entry.
    project_version = ProjectVersionDB(
        project_id=project_id,
        name="test_protection_version",
        files={file_to_keep.name: v1_details},
    )
    db_session_sync.add(project_version)
    with contextlib.suppress(IntegrityError):
        db_session_sync.commit()

    # Run the script with time mocked so deletion markers appear old enough (or not).
    mocked_now = datetime.now(timezone.utc) + timedelta(days=days_offset)
    with (
        patch("divbase_api.scripts.cleanup_expired_files.datetime") as mock_datetime_script,
        patch("divbase_api.services.s3_client.datetime") as mock_datetime_s3,
    ):
        mock_datetime_script.now.return_value = mocked_now
        mock_datetime_s3.now.return_value = mocked_now
        cleanup_expired_files.main()

    versions_purge = s3_file_manager.s3_client.list_object_versions(Bucket=bucket_name, Prefix=file_to_purge.name)
    versions_keep = s3_file_manager.s3_client.list_object_versions(Bucket=bucket_name, Prefix=file_to_keep.name)

    if expect_deleted:
        # file_to_purge should be completely gone — no versions, no markers.
        assert "Versions" not in versions_purge or not [
            v for v in versions_purge["Versions"] if v["Key"] == file_to_purge.name
        ]
        assert "DeleteMarkers" not in versions_purge or not [
            m for m in versions_purge["DeleteMarkers"] if m["Key"] == file_to_purge.name
        ]

        # file_to_keep: v1 survives, v2 is gone, and a deletion marker is the latest.
        version_ids = [v["VersionId"] for v in versions_keep.get("Versions", []) if v["Key"] == file_to_keep.name]
        assert v1_details["version_id"] in version_ids
        assert v2_details["version_id"] not in version_ids

        latest_is_marker = any(
            m["IsLatest"] for m in versions_keep.get("DeleteMarkers", []) if m["Key"] == file_to_keep.name
        )
        assert latest_is_marker is True

    else:
        # Nothing should have been touched — retention window not yet expired.
        ids = [v["VersionId"] for v in versions_keep.get("Versions", []) if v["Key"] == file_to_keep.name]
        assert v1_details["version_id"] in ids
        assert v2_details["version_id"] in ids
        assert any(v["Key"] == file_to_purge.name for v in versions_purge.get("Versions", []))
        assert any(m["Key"] == file_to_purge.name for m in versions_purge.get("DeleteMarkers", []))


def test_active_file_is_untouched(
    s3_file_manager: S3FileManager,
    project_map: dict,
    CONSTANTS: dict,
    tmp_path: Path,
):
    """
    Validate that a non-soft-deleted file is:
        - not touched.
        - no delete marker is added.
    """
    project_name = CONSTANTS["CLEANED_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    filename = f"{tmp_path.name}_active.txt"

    f = tmp_path / filename
    f.write_text("active content")
    s3_file_manager.upload_files(bucket_name=bucket_name, to_upload={filename: f})
    initial_version_id = s3_file_manager.state_of_latest_version_of_all_files(bucket_name=bucket_name)[filename][
        "version_id"
    ]

    mocked_now = datetime.now(timezone.utc) + timedelta(days=DEFAULT_SOFT_DELETED_FILES_RETENTION_DAYS + 1)
    _run_script_with_mocked_time(mocked_now)

    versions = s3_file_manager.s3_client.list_object_versions(Bucket=bucket_name, Prefix=filename)
    surviving_ids = [v["VersionId"] for v in versions.get("Versions", []) if v["Key"] == filename]
    assert initial_version_id in surviving_ids
    assert not [m for m in versions.get("DeleteMarkers", []) if m["Key"] == filename]


def test_all_versions_protected(
    db_session_sync: Session,
    s3_file_manager: S3FileManager,
    project_map: dict,
    CONSTANTS: dict,
    tmp_path: Path,
):
    """
    Validate that when every version of a file is protected:
        - none are deleted
        - latest version is a delete marker
    """
    project_name = CONSTANTS["CLEANED_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    project_id = project_map[project_name]
    filename = f"{tmp_path.name}_all_protected.txt"

    # upload + create project versions for the files.
    f = tmp_path / filename
    f.write_text("v1")
    s3_file_manager.upload_files(bucket_name=bucket_name, to_upload={filename: f})
    v1_details = s3_file_manager.state_of_latest_version_of_all_files(bucket_name=bucket_name)[filename]

    f.write_text("v2")
    s3_file_manager.upload_files(bucket_name=bucket_name, to_upload={filename: f})
    v2_details = s3_file_manager.state_of_latest_version_of_all_files(bucket_name=bucket_name)[filename]

    s3_file_manager.soft_delete_objects(objects=[filename], bucket_name=bucket_name)

    for ver_name, ver_details in [("v1", v1_details), ("v2", v2_details)]:
        pv = ProjectVersionDB(
            project_id=project_id,
            name=f"test_all_protected_{ver_name}_{tmp_path.name}",
            files={filename: ver_details},
        )
        db_session_sync.add(pv)
        db_session_sync.commit()

    mocked_now = datetime.now(timezone.utc) + timedelta(days=DEFAULT_SOFT_DELETED_FILES_RETENTION_DAYS + 1)
    _run_script_with_mocked_time(mocked_now)

    versions = s3_file_manager.s3_client.list_object_versions(Bucket=bucket_name, Prefix=filename)
    surviving_ids = [v["VersionId"] for v in versions.get("Versions", []) if v["Key"] == filename]
    assert v1_details["version_id"] in surviving_ids
    assert v2_details["version_id"] in surviving_ids
    assert any(m["IsLatest"] for m in versions.get("DeleteMarkers", []) if m["Key"] == filename)


def test_multiple_protected_versions(
    db_session_sync: Session,
    s3_file_manager: S3FileManager,
    project_map: dict,
    CONSTANTS: dict,
    tmp_path: Path,
):
    """
    Validate for a file with multiple versions:
        - only the unprotected version is deleted.
        - both protected versions survive
        - a delete marker is the latest version afterwards.
    """
    project_name = CONSTANTS["CLEANED_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    project_id = project_map[project_name]
    filename = f"{tmp_path.name}_multi_protected.txt"

    f = tmp_path / filename
    version_details = {}
    for label in ("v1", "v2", "v3"):
        f.write_text(label)
        s3_file_manager.upload_files(bucket_name=bucket_name, to_upload={filename: f})
        version_details[label] = s3_file_manager.state_of_latest_version_of_all_files(bucket_name=bucket_name)[filename]

    s3_file_manager.soft_delete_objects(objects=[filename], bucket_name=bucket_name)

    # Protect v1 and v2, leave v3 unprotected.
    for ver_name in ("v1", "v2"):
        pv = ProjectVersionDB(
            project_id=project_id,
            name=f"test_multi_{ver_name}_{tmp_path.name}",
            files={filename: version_details[ver_name]},
        )
        db_session_sync.add(pv)
        db_session_sync.commit()

    mocked_now = datetime.now(timezone.utc) + timedelta(days=DEFAULT_SOFT_DELETED_FILES_RETENTION_DAYS + 1)
    _run_script_with_mocked_time(mocked_now)

    versions = s3_file_manager.s3_client.list_object_versions(Bucket=bucket_name, Prefix=filename)
    surviving_ids = [v["VersionId"] for v in versions.get("Versions", []) if v["Key"] == filename]
    assert version_details["v1"]["version_id"] in surviving_ids
    assert version_details["v2"]["version_id"] in surviving_ids
    assert version_details["v3"]["version_id"] not in surviving_ids
    assert any(m["IsLatest"] for m in versions.get("DeleteMarkers", []) if m["Key"] == filename)


def test_results_files_hard_deleted(
    s3_file_manager: S3FileManager,
    project_map: dict,
    CONSTANTS: dict,
    tmp_path: Path,
):
    """
    Validate that query job/task files are hard deleted after the retention window:
        - all versions and delete markers removed
        - true whether the file was soft-deleted or not (both paths exercised in one test).
    """
    project_name = CONSTANTS["CLEANED_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]

    soft_deleted_filename = f"{QUERY_RESULTS_FILE_PREFIX}{tmp_path.name}_soft.vcf.gz"
    f = tmp_path / soft_deleted_filename
    f.write_text("soft-deleted results content")
    s3_file_manager.upload_files(bucket_name=bucket_name, to_upload={soft_deleted_filename: f})
    s3_file_manager.soft_delete_objects(objects=[soft_deleted_filename], bucket_name=bucket_name)

    active_filename = f"{QUERY_RESULTS_FILE_PREFIX}{tmp_path.name}_active.vcf.gz"
    f = tmp_path / active_filename
    f.write_text("active results content")
    s3_file_manager.upload_files(bucket_name=bucket_name, to_upload={active_filename: f})

    mocked_now = datetime.now(timezone.utc) + timedelta(days=DEFAULT_RESULTS_FILES_RETENTION_DAYS + 1)
    _run_script_with_mocked_time(mocked_now)

    for filename in (soft_deleted_filename, active_filename):
        versions = s3_file_manager.s3_client.list_object_versions(Bucket=bucket_name, Prefix=filename)
        assert not [v for v in versions.get("Versions", []) if v["Key"] == filename]
        assert not [m for m in versions.get("DeleteMarkers", []) if m["Key"] == filename]
