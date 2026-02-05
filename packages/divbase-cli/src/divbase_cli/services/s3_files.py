"""
Service layer for DivBase CLI S3 file operations.
"""

import sys
from pathlib import Path

import httpx
import typer

from divbase_cli.cli_exceptions import (
    FileDoesNotExistInSpecifiedVersionError,
    FilesAlreadyInProjectError,
)
from divbase_cli.services.pre_signed_urls import (
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
from divbase_lib.api_schemas.s3 import (
    FileChecksumResponse,
    ListObjectsRequest,
    ListObjectsResponse,
    ObjectDetails,
    ObjectInfoResponse,
    PreSignedDownloadResponse,
    PreSignedSinglePartUploadResponse,
    RestoreObjectsResponse,
)
from divbase_lib.divbase_constants import (
    MAX_S3_API_BATCH_SIZE,
    QUERY_RESULTS_FILE_PREFIX,
    S3_MULTIPART_UPLOAD_THRESHOLD,
)
from divbase_lib.s3_checksums import (
    MD5CheckSumFormat,
    calculate_composite_md5_s3_etag,
    calculate_md5_checksum,
    convert_checksum_hex_to_base64,
)


def list_files_command(
    divbase_base_url: str,
    project_name: str,
    prefix_filter: str | None = None,
    include_results_files: bool = False,
) -> list[ObjectDetails]:
    """
    List all files in a project optionally filtered by a prefix.
    We page through results if there are more than can be returned in a single API call.

    NOTE: The current implementation is not very efficient as we page through all results before returning any.
    Keeping simple for now as we don't expect projects to have huge numbers of files.
    But could be revisted later if performance becomes an issue.
    """
    api_route = f"v1/s3/list?project_name={project_name}"
    initial_request = ListObjectsRequest(
        prefix=prefix_filter,
        next_token=None,
    )

    response = make_authenticated_request(
        method="POST",
        divbase_base_url=divbase_base_url,
        api_route=api_route,
        json=initial_request.model_dump(),
    )
    response_data = ListObjectsResponse(**response.json())
    all_matches = response_data.objects

    # page through any remaining results
    while response_data.next_token:
        next_request = ListObjectsRequest(prefix=prefix_filter, next_token=response_data.next_token)
        response = make_authenticated_request(
            method="POST",
            divbase_base_url=divbase_base_url,
            api_route=api_route,
            json=next_request.model_dump(),
        )
        next_page_data = ListObjectsResponse(**response.json())

        all_matches.extend(next_page_data.objects)
        response_data.next_token = next_page_data.next_token

    # To enable us to both search by prefix and optionally hide/include query results files,
    # we hide DivBase query results/job files now instead of via S3.
    # so the prefix param can be used for the optional filter the user wants.
    if not include_results_files:
        all_matches = [obj for obj in all_matches if not obj.name.startswith(QUERY_RESULTS_FILE_PREFIX)]

    return all_matches


def get_file_info_command(divbase_base_url: str, project_name: str, object_name: str) -> ObjectInfoResponse:
    """Get detailed information about a specific file/object in a project."""
    response = make_authenticated_request(
        method="GET",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/s3/info?project_name={project_name}&object_name={object_name}",
    )
    return ObjectInfoResponse(**response.json())


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


def restore_objects_command(divbase_base_url: str, project_name: str, all_files: list[str]) -> RestoreObjectsResponse:
    """
    Restore soft_deleted objects in the project's bucket.
    Returns an object containing a list of the restored objects, and those that were not restored.
    """
    response = make_authenticated_request(
        method="POST",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/s3/restore?project_name={project_name}",
        json=all_files,
    )
    return RestoreObjectsResponse(**response.json())


def download_files_command(
    divbase_base_url: str,
    project_name: str,
    raw_files_input: list[str],
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

    if project_version is not None:
        offending_files = [file for file in raw_files_input if ":" in file]
        if offending_files:
            print(
                "[red] ERROR: bad Input: If you provide a global project version (using --project-version) "
                "you cannot also specify specific versions of individual files to download. \n"
                "offending files in your input: \n"
                f"{'\n'.join(offending_files)} \n"
                "Exiting..."
            )
            raise typer.Exit(1)

    if project_version:
        project_version_details = get_version_details_command(
            project_name=project_name, divbase_base_url=divbase_base_url, version_name=project_version
        )

        # check if all files specified exist for download exist at this project version
        missing_objects = [f for f in raw_files_input if f not in project_version_details.files]
        if missing_objects:
            raise FileDoesNotExistInSpecifiedVersionError(
                project_name=project_name,
                project_version=project_version,
                missing_files=missing_objects,
            )
        # we validated in CLI command that
        to_download = {file: project_version_details.files[file] for file in raw_files_input}
        json_data = [{"name": obj, "version_id": to_download[obj]} for obj in raw_files_input]
    else:
        # parse raw file inputs to see if any specific version ids were provided using format:
        # file_name:version_id (not possible when using project_version)
        json_data = []
        for file_input in raw_files_input:
            if ":" in file_input:
                name, version_id = file_input.split(sep=":", maxsplit=1)
                json_data.append({"name": name, "version_id": version_id})
            else:
                json_data.append({"name": file_input, "version_id": None})

    successful_downloads, failed_downloads = [], []
    for i in range(0, len(json_data), MAX_S3_API_BATCH_SIZE):
        batch_json_data = json_data[i : i + MAX_S3_API_BATCH_SIZE]
        response = make_authenticated_request(
            method="POST",
            divbase_base_url=divbase_base_url,
            api_route=f"v1/s3/download?project_name={project_name}",
            json=batch_json_data,
        )
        pre_signed_urls = [PreSignedDownloadResponse(**item) for item in response.json()]

        batch_download_success, batch_download_failed = download_multiple_pre_signed_urls(
            pre_signed_urls=pre_signed_urls,
            download_dir=download_dir,
            verify_checksums=verify_checksums,
        )
        successful_downloads.extend(batch_download_success)
        failed_downloads.extend(batch_download_failed)

    return DownloadOutcome(successful=successful_downloads, failed=failed_downloads)


def stream_file_command(
    divbase_base_url: str, project_name: str, file_name: str, version_id: str | None = None
) -> None:
    """Stream the contents of a single file in the project's S3 bucket to stdout."""
    json_data = [{"name": file_name, "version_id": version_id}]
    response = make_authenticated_request(
        method="POST",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/s3/download?project_name={project_name}",
        json=json_data,
    )
    pre_signed_url = PreSignedDownloadResponse(**response.json()[0]).pre_signed_url

    try:
        with httpx.stream("GET", pre_signed_url, timeout=None) as response:
            response.raise_for_status()
            for chunk in response.iter_bytes():
                sys.stdout.buffer.write(chunk)
    except BrokenPipeError:
        # This happens when the user pipes to a command that closes early
        # (e.g., `[divbase-cli stream command] | head -n 10`).
        pass


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

    files_below_threshold, files_above_threshold = [], []
    for file in all_files:
        if file.stat().st_size <= S3_MULTIPART_UPLOAD_THRESHOLD:
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
    Calculate the checksums of all local files (to be uploaded) and compares them to the checksums of all files in the project's S3 bucket.
    Raises error if any files already exist in the project's S3 bucket with identical checksums.

    Here, we are catching an attempt to upload an identical file twice.
    This is only ran if 'safe_mode' is enabled for uploads.
    We do not catch an attempt to upload an identical object if it has a different name.

    Return a dict of file names with hex-encoded checksums for all files to be uploaded (including those that are not in S3).
    These checksums are later used when uploading to the server so the server can verify the upload.
    """
    already_uploaded_files: dict[Path, str] = {}  # files that already exist in S3 with identical checksum
    local_checksums: dict[str, str] = {}  # all local files checksums

    # have to batch requests if above max number allowed by divbase server
    for i in range(0, len(all_files), MAX_S3_API_BATCH_SIZE):
        batch_files = all_files[i : i + MAX_S3_API_BATCH_SIZE]
        batch_files_names = [file.name for file in batch_files]

        response = make_authenticated_request(
            method="POST",
            divbase_base_url=divbase_base_url,
            api_route=f"v1/s3/checksums?project_name={project_name}",
            json=batch_files_names,
        )
        server_checksum_responses = [FileChecksumResponse(**item) for item in response.json()]
        server_checksums = {item.object_name: item.md5_checksum for item in server_checksum_responses}

        for file in batch_files:
            if file.stat().st_size > S3_MULTIPART_UPLOAD_THRESHOLD:
                calculated_checksum = calculate_composite_md5_s3_etag(file_path=file)
            else:
                calculated_checksum = calculate_md5_checksum(file_path=file, output_format=MD5CheckSumFormat.HEX)

            local_checksums[file.name] = calculated_checksum
            if server_checksums.get(file.name) and server_checksums[file.name] == calculated_checksum:
                already_uploaded_files[file] = calculated_checksum

    if already_uploaded_files:
        raise FilesAlreadyInProjectError(existing_files=already_uploaded_files, project_name=project_name)

    return local_checksums
