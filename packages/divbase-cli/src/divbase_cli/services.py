"""
Service layer for DivBase CLI project version and S3 file operations.
"""

from pathlib import Path

from divbase_cli.cli_exceptions import (
    FileDoesNotExistInSpecifiedVersionError,
    FilesAlreadyInProjectError,
)
from divbase_cli.pre_signed_urls import (
    DownloadOutcome,
    UploadOutcome,
    download_multiple_pre_signed_urls,
    upload_multiple_pre_signed_urls,
)
from divbase_cli.user_auth import make_authenticated_request
from divbase_lib.api_schemas.project_versions import (
    AddVersionRequest,
    AddVersionResponse,
    DeleteVersionRequest,
    DeleteVersionResponse,
    ProjectVersionDetailResponse,
    ProjectVersionInfo,
)
from divbase_lib.api_schemas.s3 import ExistingFileResponse, PreSignedDownloadResponse, PreSignedUploadResponse
from divbase_lib.s3_checksums import MD5CheckSumFormat, calculate_md5_checksum, convert_checksum_hex_to_base64


def add_version_command(project_name: str, divbase_base_url: str, name: str, description: str) -> AddVersionResponse:
    """Add a new version to the project versions table stored on the divbase server"""
    request_data = AddVersionRequest(name=name, description=description)

    response = make_authenticated_request(
        method="PATCH",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/project-versions/add?project_name={project_name}",
        json=request_data.model_dump(),
    )

    return AddVersionResponse(**response.json())


def list_versions_command(project_name: str, include_deleted: bool, divbase_base_url: str) -> list[ProjectVersionInfo]:
    """
    List all versions in the project versions table stored on the divbase server.
    Returns a dict of version names (keys) to details about the versions.
    """
    response = make_authenticated_request(
        method="GET",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/project-versions/list?project_name={project_name}&include_deleted={str(include_deleted).lower()}",
    )

    project_versions = []
    response_data = response.json()
    for version in response_data:
        project_versions.append(ProjectVersionInfo(**version))

    return project_versions


def get_version_details_command(
    project_name: str, divbase_base_url: str, version_name: str
) -> ProjectVersionDetailResponse:
    """Get details about a specific project version, including all files and their version IDs at that version."""
    response = make_authenticated_request(
        method="GET",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/project-versions/version_details?project_name={project_name}&version_name={version_name}",
    )

    return ProjectVersionDetailResponse(**response.json())


def delete_version_command(project_name: str, divbase_base_url: str, version_name: str) -> DeleteVersionResponse:
    """
    Delete a version from the project versions table stored on the divbase server.
    This marks the version as (soft) deleted server side,
    and it will eventually be permanently deleted (after some grace period).
    """
    request_data = DeleteVersionRequest(version_name=version_name)

    response = make_authenticated_request(
        method="DELETE",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/project-versions/delete?project_name={project_name}",
        json=request_data.model_dump(),
    )

    return DeleteVersionResponse(**response.json())


def list_files_command(divbase_base_url: str, project_name: str) -> list[str]:
    """List all files in a project."""
    response = make_authenticated_request(
        method="GET",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/s3/?project_name={project_name}",
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
    file_checksums_hex = {}
    if safe_mode:
        for file in all_files:
            file_checksums_hex[file.name] = calculate_md5_checksum(file_path=file, output_format=MD5CheckSumFormat.HEX)

        files_to_check = []
        for file in all_files:
            files_to_check.append({"object_name": file.name, "md5_checksum": file_checksums_hex[file.name]})

        response = make_authenticated_request(
            method="POST",
            divbase_base_url=divbase_base_url,
            api_route=f"v1/s3/check-exists?project_name={project_name}",
            json=files_to_check,
        )
        existing_files = response.json()

        if existing_files:
            existing_object_names = [ExistingFileResponse(**file) for file in existing_files]
            raise FilesAlreadyInProjectError(existing_files=existing_object_names, project_name=project_name)

    objects_to_upload = []
    for file in all_files:
        upload_object = {
            "name": file.name,
            "content_length": file.stat().st_size,
        }
        if safe_mode:
            hex_checksum = file_checksums_hex[file.name]
            base64_checksum = convert_checksum_hex_to_base64(hex_checksum)
            upload_object["md5_hash"] = base64_checksum

        objects_to_upload.append(upload_object)

    response = make_authenticated_request(
        method="POST",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/s3/upload?project_name={project_name}",
        json=objects_to_upload,
    )
    pre_signed_urls = [PreSignedUploadResponse(**item) for item in response.json()]
    return upload_multiple_pre_signed_urls(pre_signed_urls=pre_signed_urls, all_files=all_files)


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
