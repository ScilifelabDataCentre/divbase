"""
Service layer for DivBase CLI S3 file operations.
"""

import sys
from dataclasses import dataclass
from pathlib import Path

import httpx
import typer
from rich import print

from divbase_cli.cli_exceptions import (
    FileDoesNotExistInSpecifiedVersionError,
)
from divbase_cli.services.pre_signed_urls import (
    DownloadOutcome,
    FailedUpload,
    SkippedUpload,
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
    MakeDirectoriesResponse,
    ObjectDetails,
    ObjectInfoResponse,
    PreSignedDownloadResponse,
    PreSignedSinglePartUploadResponse,
    RestoreObjectsResponse,
    SoftDeletedObjectDetails,
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


@dataclass
class ToDownload:
    """
    Represent a file to be download in the download-all command. This unifies the 2 ways we get the files needed:
    1. From the latest files in the bucket (listing files)
    2. From a user defined project version (getting the files from the version details)"
    """

    name: str
    etag: str
    size_bytes: int
    version_id: str | None = None  # latest version if None


@dataclass
class ToUpload:
    """
    Represents a file to be uploaded in the files upload command.
    """

    file_path: Path
    destination_key: str  # what the file will be called in the project store. Can include "/"s to represent "folder" paths e.g. "some/subdirs/file.tsv"
    checksum_local: str | None = None  # None if not uploaded with "safe-mode"

    @property
    def file_name(self) -> str:
        return self.file_path.name

    @property
    def file_size(self) -> int:
        return self.file_path.stat().st_size


def list_files_command(
    divbase_base_url: str,
    project_name: str,
    prefix: str | None = None,
    include_results_files: bool = False,
    file_system_view: bool = True,
) -> tuple[list[ObjectDetails], list[str]]:
    """
    List files and (optionally) folders in a project.

    file_system_view = true will return both files and folders at the current prefix level, simulating a file system view.
    Otherwise all files are returned in a flat list
    TODO - consider that we page through all results before returning.
    """
    if file_system_view:
        delimiter = "/"
    else:
        delimiter = None

    api_route = f"v1/s3/list?project_name={project_name}"
    initial_request = ListObjectsRequest(prefix=prefix, delimiter=delimiter, next_token=None)

    response = make_authenticated_request(
        method="POST",
        divbase_base_url=divbase_base_url,
        api_route=api_route,
        json=initial_request.model_dump(),
    )
    response_data = ListObjectsResponse(**response.json())
    folders, files = list(response_data.folders), list(response_data.files)

    while response_data.next_token:
        next_request = ListObjectsRequest(prefix=prefix, delimiter=delimiter, next_token=response_data.next_token)
        response = make_authenticated_request(
            method="POST",
            divbase_base_url=divbase_base_url,
            api_route=api_route,
            json=next_request.model_dump(),
        )
        response_data = ListObjectsResponse(**response.json())
        folders.extend(response_data.folders)
        files.extend(response_data.files)

    # Hide DivBase query results files unless explicitly requested.
    # we do this here instead of in the api query so we can use the prefix filter to filter on other things
    if not include_results_files:
        files = [obj for obj in files if not obj.name.startswith(QUERY_RESULTS_FILE_PREFIX)]

    return files, folders


def list_soft_deleted_files_command(divbase_base_url: str, project_name: str) -> list[SoftDeletedObjectDetails]:
    """List all soft-deleted files in a project."""
    api_route = f"v1/s3/list/soft-deleted?project_name={project_name}"
    response = make_authenticated_request(
        method="GET",
        divbase_base_url=divbase_base_url,
        api_route=api_route,
    )
    return [SoftDeletedObjectDetails(**obj) for obj in response.json()]


def get_file_info_command(divbase_base_url: str, project_name: str, object_name: str) -> ObjectInfoResponse:
    """Get detailed information about a specific file/object in a project."""
    response = make_authenticated_request(
        method="GET",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/s3/info?project_name={project_name}&object_name={object_name}",
    )
    return ObjectInfoResponse(**response.json())


def make_directories_command(
    divbase_base_url: str, project_name: str, directories: list[str]
) -> MakeDirectoriesResponse:
    """
    Create directories in the project's store.
    Pagination not supported here, so up to MAX_S3_API_BATCH_SIZE directories can be created at once.
    """
    response = make_authenticated_request(
        method="POST",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/s3/mkdir?project_name={project_name}",
        json=directories,
    )
    return MakeDirectoriesResponse(**response.json())


def soft_delete_objects_command(divbase_base_url: str, project_name: str, all_files: list[str]) -> list[str]:
    """
    Soft delete objects from the project's bucket.
    Returns a list of the soft deleted objects
    """
    deleted_objects = []
    for i in range(0, len(all_files), MAX_S3_API_BATCH_SIZE):
        batch_files = all_files[i : i + MAX_S3_API_BATCH_SIZE]
        response = make_authenticated_request(
            method="DELETE",
            divbase_base_url=divbase_base_url,
            api_route=f"v1/s3/?project_name={project_name}",
            json=batch_files,
        )
        deleted_objects.extend(response.json())
    return deleted_objects


def restore_objects_command(divbase_base_url: str, project_name: str, all_files: list[str]) -> RestoreObjectsResponse:
    """
    Restore soft_deleted objects in the project's bucket.
    Returns an object containing a list of the restored objects, and those that were not restored.
    """
    all_restored = []
    all_not_restored = []
    for i in range(0, len(all_files), MAX_S3_API_BATCH_SIZE):
        batch_files = all_files[i : i + MAX_S3_API_BATCH_SIZE]
        response = make_authenticated_request(
            method="POST",
            divbase_base_url=divbase_base_url,
            api_route=f"v1/s3/restore?project_name={project_name}",
            json=batch_files,
        )
        batch_response = RestoreObjectsResponse(**response.json())
        all_restored.extend(batch_response.restored)
        all_not_restored.extend(batch_response.not_restored)

    return RestoreObjectsResponse(restored=all_restored, not_restored=all_not_restored)


def download_files_command(
    divbase_base_url: str,
    project_name: str,
    raw_files_input: list[str],
    download_dir: Path,
    verify_checksums: bool,
    dry_run: bool,
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

        json_data = []
        for file_name, file_details in project_version_details.files.items():
            if file_name in raw_files_input:
                json_data.append({"name": file_name, "version_id": file_details["version_id"]})
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

    if dry_run:
        print("\n[green bold]Dry run enabled, The following files would have been downloaded:[/green bold]")
        for file_info in json_data:
            version_id = file_info["version_id"]
            if version_id:
                print(f"- '{file_info['name']}' (version_id: '{version_id}')")
            else:
                print(f"- '{file_info['name']}' (latest version)")
        raise typer.Exit(0)

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
    project_name: str,
    divbase_base_url: str,
    all_files: list[ToUpload],
    safe_mode: bool,
    skip_existing: bool = False,
    dry_run: bool = False,
) -> UploadOutcome:
    """
    Upload files to the project's S3 bucket.
    Returns an UploadOutcome object containing details about how uploading went for each file.

    - Safe mode:
        1. checks if any of the files that are to be uploaded already exist in the bucket (by comparing checksums)
        2. Adds checksum to upload request to allow server to verify upload.
    """
    all_successful_uploads: list[SuccessfulUpload] = []
    all_failed_uploads: list[FailedUpload] = []
    all_skipped_uploads: list[SkippedUpload] = []

    if safe_mode:
        if dry_run:
            print(
                "[green]Dry run enabled, calculating the checksums of all files to upload to determine which would be uploaded vs skipped...[/green]"
            )
        else:
            print("[green]Preparing to upload files, calculating the checksums of all files to upload...[/green]")

        to_upload, already_uploaded = filter_already_uploaded_files(
            project_name=project_name,
            divbase_base_url=divbase_base_url,
            all_files=all_files,
        )

        if already_uploaded:
            # skip_existing means we just skip uploading those files that are already there
            if not skip_existing:
                if dry_run:
                    print(
                        "[red bold]Error: The following upload attempt would have failed due to the below error:\n[/red bold]"
                    )

                files_str = "\n".join(f"'{f.file_path}' (Checksum: {f.checksum_local})" for f in already_uploaded)
                print(
                    f"\n[red bold]Error: For the project: '{project_name}'\n"
                    "The exact version of the following files that you're trying to upload already exist inside the project:\n[/red bold]"
                    f"[red]{files_str}[/red]\n"
                    "[bold green]Tip: if you want to skip re-uploading these files and continue uploading the other files, re-run this command with the '--skip-existing' flag.[/bold green]"
                )
                raise typer.Exit(1)
            else:
                for file in already_uploaded:
                    all_skipped_uploads.append(
                        SkippedUpload(
                            object_name=file.destination_key,
                            file_path=file.file_path,
                            reason=f"The file with checksum {file.checksum_local} already exists in the project",
                        )
                    )
    else:
        # we don't provide checksums to server
        to_upload = all_files

    if dry_run:
        if to_upload:
            print("\n[green bold]The following files would have been uploaded:[/green bold]")
            for file in to_upload:
                print(f"- '{file.file_path}' -> would be stored as: '{file.destination_key}' in the project")
        if all_skipped_uploads:
            print("\n[yellow bold]The following files would have been skipped from upload:[/yellow bold]")
            for file in all_skipped_uploads:
                print(f"- [yellow]'{file.file_path}' Reason: {file.reason}[/yellow]")

        raise typer.Exit(0)

    # Split files into those that need single vs multipart upload
    files_below_threshold: list[ToUpload] = []
    files_above_threshold: list[ToUpload] = []
    for file in to_upload:
        if file.file_size <= S3_MULTIPART_UPLOAD_THRESHOLD:
            files_below_threshold.append(file)
        else:
            files_above_threshold.append(file)

    # P1. Process all single-part uploads in batches of max size allowed by divbase server.
    for i in range(0, len(files_below_threshold), MAX_S3_API_BATCH_SIZE):
        batch_files = files_below_threshold[i : i + MAX_S3_API_BATCH_SIZE]
        batch_of_objects_to_upload = []
        for file in batch_files:
            upload_object = {
                "name": file.destination_key,
                "content_length": file.file_size,
            }
            if safe_mode and file.checksum_local:
                # server expects base64 encoded checksum when provided
                upload_object["md5_hash"] = convert_checksum_hex_to_base64(hex_checksum=file.checksum_local)
            batch_of_objects_to_upload.append(upload_object)

        response = make_authenticated_request(
            method="POST",
            divbase_base_url=divbase_base_url,
            api_route=f"v1/s3/upload/single-part?project_name={project_name}",
            json=batch_of_objects_to_upload,
        )
        pre_signed_urls = [PreSignedSinglePartUploadResponse(**item) for item in response.json()]
        dest_key_to_path = {f.destination_key: f.file_path for f in batch_files}

        successful_uploads, failed_uploads = upload_multiple_singlepart_pre_signed_urls(
            pre_signed_urls=pre_signed_urls, dest_key_to_path=dest_key_to_path
        )
        all_successful_uploads.extend(successful_uploads)
        all_failed_uploads.extend(failed_uploads)

    # P2. process all multipart uploads.
    for file in files_above_threshold:
        outcome = perform_multipart_upload(
            project_name=project_name,
            divbase_base_url=divbase_base_url,
            file_path=file.file_path,
            safe_mode=safe_mode,
            destination_name=file.destination_key,
        )

        if isinstance(outcome, SuccessfulUpload):
            all_successful_uploads.append(outcome)
        else:
            all_failed_uploads.append(outcome)

    return UploadOutcome(successful=all_successful_uploads, failed=all_failed_uploads, skipped=all_skipped_uploads)


def filter_already_uploaded_files(
    project_name: str,
    divbase_base_url: str,
    all_files: list[ToUpload],
) -> tuple[list[ToUpload], list[ToUpload]]:
    """
    Separate files into whether they are already in the projects bucket or not.
    (1st list to be uploaded, 2nd list is files already in projects bucket)

    For a file to be considered already in the project,
    there must be an object with the same key and checksum in the project's S3 bucket.

    This is only ran if 'safe_mode' is enabled for uploads.
    These checksums are later used when uploading to the server so the server can verify the upload.
    """
    already_uploaded: list[ToUpload] = []
    to_be_uploaded: list[ToUpload] = []

    # have to batch requests if above max number allowed by divbase server
    for i in range(0, len(all_files), MAX_S3_API_BATCH_SIZE):
        batch_files = all_files[i : i + MAX_S3_API_BATCH_SIZE]
        batch_object_keys = [file.destination_key for file in batch_files]

        response = make_authenticated_request(
            method="POST",
            divbase_base_url=divbase_base_url,
            api_route=f"v1/s3/checksums?project_name={project_name}",
            json=batch_object_keys,
        )
        server_checksum_responses = [FileChecksumResponse(**item) for item in response.json()]
        server_checksums = {item.object_name: item.md5_checksum for item in server_checksum_responses}

        for file in batch_files:
            local_checksum = _calc_local_checksum(file_path=file.file_path)
            file_to_upload = ToUpload(
                file_path=file.file_path, destination_key=file.destination_key, checksum_local=local_checksum
            )

            if server_checksums.get(file.destination_key) and server_checksums[file.destination_key] == local_checksum:
                already_uploaded.append(file_to_upload)
            else:
                to_be_uploaded.append(file_to_upload)

            print(f"MD5 Checksum calculated for file: '{file.destination_key}'")

    return to_be_uploaded, already_uploaded


def filter_out_already_downloaded_files(
    all_files: list[ToDownload], download_dir: Path
) -> tuple[list[ToDownload], list[ToDownload]]:
    """
    Filter out files that already exist in a local directory with the same checksum.

    Returns two lists:
    1. Files that do not exist locally and need to be downloaded.
    2. Files that are not identical to the file in S3 and will be overwritten if the user wants to download them.
    """
    files_to_download, files_to_overwrite = [], []

    for s3_file in all_files:
        local_file_path = download_dir / s3_file.name

        if not local_file_path.exists():
            files_to_download.append(s3_file)
            continue

        local_checksum = _calc_local_checksum(file_path=local_file_path)
        if local_checksum != s3_file.etag:
            files_to_overwrite.append(s3_file)

    return files_to_download, files_to_overwrite


def _calc_local_checksum(file_path: Path) -> str:
    """
    Calculate the checksum for a local file. Handles whether to use a single or composite checksum based on file size.

    This is used to validate the integrity of a file before upload.
    Or determine if a file can be skipped from downloading.
    """
    if file_path.stat().st_size > S3_MULTIPART_UPLOAD_THRESHOLD:
        return calculate_composite_md5_s3_etag(file_path=file_path)
    else:
        return calculate_md5_checksum(file_path=file_path, output_format=MD5CheckSumFormat.HEX)
