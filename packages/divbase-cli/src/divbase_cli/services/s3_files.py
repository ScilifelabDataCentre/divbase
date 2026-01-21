"""
Service layer for DivBase CLI S3 file operations.
"""

from pathlib import Path

from divbase_cli.cli_exceptions import (
    FileDoesNotExistInSpecifiedVersionError,
    FilesAlreadyInProjectError,
)
from divbase_cli.services.pre_signed_urls import (
    MULTIPART_UPLOAD_THRESHOLD,
    DownloadOutcome,
    FailedUpload,
    SuccessfulUpload,
    UploadOutcome,
    download_multiple_pre_signed_urls,
    perform_multipart_upload,
    upload_multiple_singlepart_pre_signed_urls,
)
from divbase_cli.services.project_versions import get_version_details_command
from divbase_cli.user_auth import make_authenticated_request
from divbase_lib.api_schemas.divbase_constants import MAX_S3_API_BATCH_SIZE
from divbase_lib.api_schemas.s3 import (
    ExistingFileResponse,
    PreSignedDownloadResponse,
    PreSignedSinglePartUploadResponse,
)
from divbase_lib.s3_checksums import MD5CheckSumFormat, calculate_md5_checksum, convert_checksum_hex_to_base64


def list_files_command(divbase_base_url: str, project_name: str) -> list[str]:
    """List all files in a project."""
    response = make_authenticated_request(
        method="GET",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/s3/?project_name={project_name}",
    )

    return response.json()


def soft_delete_objects_command(divbase_base_url: str, project_name: str, all_files: list[str]) -> list[str]:
    """
    Soft delete objects from the project's bucket.
    Returns a list of the soft deleted objects
    """
    response = make_authenticated_request(
        method="DELETE",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/s3/?project_name={project_name}",
        json=all_files,
    )
    return response.json()


def download_files_command(
    divbase_base_url: str,
    project_name: str,
    all_files: list[str],
    download_dir: Path,
    verify_checksums: bool,
    project_version: str | None = None,
) -> DownloadOutcome:
    """
    Download files from the given project's S3 bucket.
    """
    if not download_dir.is_dir():
        raise NotADirectoryError(
            f"The specified download directory '{download_dir}' is not a directory. Please create it or specify a valid directory before continuing."
        )

    if project_version:
        project_version_details = get_version_details_command(
            project_name=project_name, divbase_base_url=divbase_base_url, version_name=project_version
        )

        # check if all files specified exist for download exist at this project version
        missing_objects = [f for f in all_files if f not in project_version_details.files]
        if missing_objects:
            raise FileDoesNotExistInSpecifiedVersionError(
                project_name=project_name,
                project_version=project_version,
                missing_files=missing_objects,
            )
        to_download = {file: project_version_details.files[file] for file in all_files}
        json_data = [{"name": obj, "version_id": to_download[obj]} for obj in all_files]
    else:
        json_data = [{"name": obj, "version_id": None} for obj in all_files]

    response = make_authenticated_request(
        method="POST",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/s3/download?project_name={project_name}",
        json=json_data,
    )
    pre_signed_urls = [PreSignedDownloadResponse(**item) for item in response.json()]

    download_results = download_multiple_pre_signed_urls(
        pre_signed_urls=pre_signed_urls, download_dir=download_dir, verify_checksums=verify_checksums
    )
    return download_results


def upload_files_command(
    project_name: str, divbase_base_url: str, all_files: list[Path], safe_mode: bool
) -> UploadOutcome:
    """
    Upload files to the project's S3 bucket.
    Returns an UploadOutcome object containing details of which files were successfully uploaded and which failed.

    - Safe mode:
        1. checks if any of the files that are to be uploaded already exist in the bucket (by comparing checksums)
        2. Adds checksum to upload request to allow server to verify upload.
    """
    if safe_mode:
        # mapping of file name to hex-encoded checksum
        file_checksums_hex = compare_local_to_s3_checksums(
            project_name=project_name,
            divbase_base_url=divbase_base_url,
            all_files=all_files,
        )
    else:
        file_checksums_hex = {}

    files_below_threshold, files_above_threshold = [], []
    for file in all_files:
        if file.stat().st_size <= MULTIPART_UPLOAD_THRESHOLD:
            files_below_threshold.append(file)
        else:
            files_above_threshold.append(file)

    all_successful_uploads: list[SuccessfulUpload] = []
    all_failed_uploads: list[FailedUpload] = []

    # P1. Process all single-part uploads in batches of max size allowed by divbase server.
    for i in range(0, len(files_below_threshold), MAX_S3_API_BATCH_SIZE):
        batch_files = files_below_threshold[i : i + MAX_S3_API_BATCH_SIZE]
        batch_of_objects_to_upload = []
        for file in batch_files:
            upload_object = {
                "name": file.name,
                "content_length": file.stat().st_size,
            }
            if safe_mode:
                hex_checksum = file_checksums_hex[file.name]
                upload_object["md5_hash"] = convert_checksum_hex_to_base64(hex_checksum)
            batch_of_objects_to_upload.append(upload_object)

        response = make_authenticated_request(
            method="POST",
            divbase_base_url=divbase_base_url,
            api_route=f"v1/s3/upload/single-part?project_name={project_name}",
            json=batch_of_objects_to_upload,
        )
        pre_signed_urls = [PreSignedSinglePartUploadResponse(**item) for item in response.json()]
        single_part_upload_outcome = upload_multiple_singlepart_pre_signed_urls(
            pre_signed_urls=pre_signed_urls, all_files=batch_files
        )
        all_successful_uploads.extend(single_part_upload_outcome.successful)
        all_failed_uploads.extend(single_part_upload_outcome.failed)

    # P2. process all multipart uploads.
    for file_path in files_above_threshold:
        outcome = perform_multipart_upload(
            project_name=project_name,
            divbase_base_url=divbase_base_url,
            file_path=file_path,
            safe_mode=safe_mode,
        )

        if isinstance(outcome, SuccessfulUpload):
            all_successful_uploads.append(outcome)
        else:
            all_failed_uploads.append(outcome)

    return UploadOutcome(successful=all_successful_uploads, failed=all_failed_uploads)


def compare_local_to_s3_checksums(project_name: str, divbase_base_url: str, all_files: list[Path]) -> dict[str, str]:
    """
    Calculate the checksums of all local files (to be uploaded) and compare them to the files (with same name - if they exist)
    already in the project's S3 bucket.

    Here we are catching an attempt to upload the same file with the exact same content (checksum) twice.
    Only ran if 'safe_mode' is enabled for uploads.

    Return a dict of file names to hex-encoded checksums for all files to be uploaded.
    """
    file_checksums_hex = {}
    for file in all_files:
        file_checksums_hex[file.name] = calculate_md5_checksum(file_path=file, output_format=MD5CheckSumFormat.HEX)

    files_to_check = []
    for file in all_files:
        files_to_check.append({"object_name": file.name, "md5_checksum": file_checksums_hex[file.name]})

    # api accepts up to 100 files to check at a time
    existing_files = []
    for i in range(0, len(files_to_check), MAX_S3_API_BATCH_SIZE):
        batch = files_to_check[i : i + MAX_S3_API_BATCH_SIZE]
        response = make_authenticated_request(
            method="POST",
            divbase_base_url=divbase_base_url,
            api_route=f"v1/s3/check-exists?project_name={project_name}",
            json=batch,
        )
        existing_files.extend(response.json())

    if existing_files:
        existing_object_names = [ExistingFileResponse(**file) for file in existing_files]
        raise FilesAlreadyInProjectError(existing_files=existing_object_names, project_name=project_name)

    return file_checksums_hex
