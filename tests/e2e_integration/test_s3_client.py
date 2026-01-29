from pathlib import Path

import pytest

from divbase_api.services.s3_client import S3FileManager, create_s3_file_manager
from divbase_lib.divbase_constants import S3_MULTIPART_UPLOAD_THRESHOLD
from divbase_lib.exceptions import ObjectDoesNotExistError


@pytest.fixture(scope="module")
def s3_client(CONSTANTS) -> S3FileManager:
    """Provides an S3FileManager instance configured for the test environment."""
    return create_s3_file_manager(url=CONSTANTS["MINIO_URL"])


@pytest.fixture(scope="module")
def large_file(tmp_path_factory) -> Path:
    """Creates a dummy file large enough to trigger a multipart upload."""
    file_path = tmp_path_factory.mktemp("large_files") / "large_file.bin"
    # Create a file slightly larger than the multipart threshold
    file_size = S3_MULTIPART_UPLOAD_THRESHOLD + 1024
    with open(file_path, "wb") as f:
        f.write(b"\0" * file_size)
    return file_path


def test_upload_and_download_one_file_singlepart(s3_client: S3FileManager, tmp_path: Path, CONSTANTS):
    """
    Tests uploading and then downloading a single, small file.
    """
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"]["project1"]
    file_to_upload = Path("tests/fixtures/file1.txt")
    object_name = file_to_upload.name
    download_dir = tmp_path / "downloads"
    download_dir.mkdir(parents=True, exist_ok=True)

    s3_client.upload_files(
        to_upload={object_name: file_to_upload},
        bucket_name=bucket_name,
    )

    s3_client.download_files(
        bucket_name=bucket_name,
        objects={object_name: None},
        download_dir=download_dir,
    )
    assert (download_dir / object_name).exists()
    assert (download_dir / object_name).read_text() == file_to_upload.read_text()


def test_upload_and_download_multiple_files_singlepart(s3_client: S3FileManager, tmp_path: Path, CONSTANTS):
    """
    Tests uploading and then downloading multiple small files.
    """
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"]["project1"]
    files_to_upload = {
        "file1.txt": Path("tests/fixtures/file1.txt"),
        "file2.txt": Path("tests/fixtures/file2.txt"),
        "file3.txt": Path("tests/fixtures/file3.txt"),
    }
    download_dir = tmp_path / "downloads"
    download_dir.mkdir(parents=True, exist_ok=True)

    s3_client.upload_files(to_upload=files_to_upload, bucket_name=bucket_name)

    s3_client.download_files(
        bucket_name=bucket_name,
        objects={name: None for name in files_to_upload},
        download_dir=download_dir,
    )
    for object_name, file_to_upload in files_to_upload.items():
        assert (download_dir / object_name).exists()
        assert (download_dir / object_name).read_text() == file_to_upload.read_text()


def test_upload_and_download_multipart(s3_client: S3FileManager, tmp_path: Path, large_file: Path, CONSTANTS):
    """
    Tests uploading and then downloading a large file, triggering multipart logic.
    """
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"]["project1"]
    object_name = large_file.name

    s3_client.upload_files(
        to_upload={object_name: large_file},
        bucket_name=bucket_name,
    )

    s3_client.download_files(
        bucket_name=bucket_name,
        objects={object_name: None},
        download_dir=tmp_path,
    )
    downloaded_path = tmp_path / object_name
    assert downloaded_path.exists()
    assert downloaded_path.stat().st_size == large_file.stat().st_size


def test_download_s3_file_to_str(s3_client: S3FileManager, CONSTANTS):
    """
    Tests downloading a file's content directly to a string.
    """
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"]["project1"]
    file_to_upload = Path("tests/fixtures/file1.txt")
    object_name = "download_to_str.txt"

    s3_client.upload_files(to_upload={object_name: file_to_upload}, bucket_name=bucket_name)

    content = s3_client.download_s3_file_to_str(key=object_name, bucket_name=bucket_name)
    assert content == file_to_upload.read_text()


def test_download_non_existent_file(s3_client: S3FileManager, tmp_path: Path, CONSTANTS):
    """
    Tests that downloading a non-existent file raises ObjectDoesNotExistError.
    """
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"]["project1"]
    object_name = "this_file_does_not_exist.txt"

    with pytest.raises(ObjectDoesNotExistError):
        s3_client.download_files(
            bucket_name=bucket_name,
            objects={object_name: None},
            download_dir=tmp_path,
        )


def test_list_files(s3_client: S3FileManager, CONSTANTS):
    """
    Tests listing files in a bucket.
    """
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"]["project1"]
    files_to_upload = {
        "list_test_1.txt": Path("tests/fixtures/file1.txt"),
        "list_test_2.txt": Path("tests/fixtures/file2.txt"),
    }

    s3_client.upload_files(to_upload=files_to_upload, bucket_name=bucket_name)

    listed_files = s3_client.list_files(bucket_name=bucket_name)
    for object_name in files_to_upload:
        assert object_name in listed_files


def test_soft_delete_objects(s3_client: S3FileManager, CONSTANTS, tmp_path: Path):
    """
    Tests that soft-deleting objects works correctly.
    """
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"]["project1"]
    files_to_upload = {
        "file_to_delete_1.txt": Path("tests/fixtures/file1.txt"),
        "file_to_delete_2.txt": Path("tests/fixtures/file2.txt"),
    }
    object_names = list(files_to_upload.keys())

    s3_client.upload_files(to_upload=files_to_upload, bucket_name=bucket_name)

    deleted_objects = s3_client.soft_delete_objects(objects=object_names, bucket_name=bucket_name)
    assert sorted(deleted_objects) == sorted(object_names)

    with pytest.raises(ObjectDoesNotExistError):
        s3_client.download_files(
            bucket_name=bucket_name,
            objects={object_names[0]: None},
            download_dir=tmp_path,
        )


def test_versioning_and_hard_delete(s3_client: S3FileManager, CONSTANTS):
    """
    Tests file versioning by uploading twice, getting the latest version,
    and then performing a hard delete on a specific version.
    """
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"]["project1"]
    object_name = "versioned_file.txt"
    file_v1_path = Path("tests/fixtures/file1.txt")
    file_v2_path = Path("tests/fixtures/file2.txt")

    s3_client.upload_files(to_upload={object_name: file_v1_path}, bucket_name=bucket_name)
    versions_after_v1 = s3_client.latest_version_of_all_files(bucket_name=bucket_name)
    v1_id = versions_after_v1[object_name]

    s3_client.upload_files(to_upload={object_name: file_v2_path}, bucket_name=bucket_name)
    versions_after_v2 = s3_client.latest_version_of_all_files(bucket_name=bucket_name)
    v2_id = versions_after_v2[object_name]

    assert v1_id != v2_id

    # Hard delete the original version and verify it no longer exists
    s3_client.hard_delete_specific_object_versions(versioned_objects={object_name: v1_id}, bucket_name=bucket_name)
    with pytest.raises(ObjectDoesNotExistError):
        s3_client.download_files(
            bucket_name=bucket_name,
            objects={object_name: v1_id},
            download_dir=Path("/tmp"),
        )

    # Verify the latest version still has the id of the second upload
    final_versions = s3_client.latest_version_of_all_files(bucket_name=bucket_name)
    assert final_versions[object_name] == v2_id
