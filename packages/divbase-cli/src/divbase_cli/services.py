"""
CLI commands for managing S3 bucket versions.
"""

from pathlib import Path
from urllib.parse import urlencode

from divbase_cli.pre_signed_urls import download_multiple_pre_signed_urls, upload_multiple_pre_signed_urls
from divbase_cli.user_auth import make_authenticated_request
from divbase_cli.user_config import ProjectConfig
from divbase_lib.exceptions import FilesAlreadyInBucketError
from divbase_lib.s3_client import create_s3_file_manager
from divbase_lib.schemas.bucket_versions import (
    AddVersionRequest,
    AddVersionResponse,
    CreateVersioningFileRequest,
    CreateVersioningFileResponse,
    DeleteVersionRequest,
    DeleteVersionResponse,
    FilesAtVersionResponse,
    VersionListResponse,
)
from divbase_lib.vcf_dimension_indexing import VCFDimensionIndexManager


def create_version_object_command(
    project_name: str, divbase_base_url: str, version_name: str = "v0.0.0", description: str = "First version"
) -> CreateVersioningFileResponse:
    """Create the initial bucket versioning file for a project."""
    request_data = CreateVersioningFileRequest(name=version_name, description=description)

    response = make_authenticated_request(
        method="POST",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/bucket-versions/create?project_name={project_name}",
        json_data=request_data.model_dump(),
    )

    return CreateVersioningFileResponse(**response.json())


def add_version_command(
    project_name: str, divbase_base_url: str, name: str, description: str = ""
) -> AddVersionResponse:
    """Add a new version to the bucket versioning file"""
    request_data = AddVersionRequest(name=name, description=description)

    response = make_authenticated_request(
        method="PATCH",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/bucket-versions/add?project_name={project_name}",
        json_data=request_data.model_dump(),
    )

    return AddVersionResponse(**response.json())


def list_versions_command(project_name: str, divbase_base_url: str) -> VersionListResponse:
    """List all versions in the bucket versioning file"""
    response = make_authenticated_request(
        method="GET",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/bucket-versions/list?project_name={project_name}",
    )

    return VersionListResponse(**response.json())


def list_files_at_version_command(
    project_name: str, divbase_base_url: str, bucket_version: str
) -> FilesAtVersionResponse:
    """List all files at a specific version"""
    response = make_authenticated_request(
        method="GET",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/bucket-versions/list_detailed?project_name={project_name}&bucket_version={bucket_version}",
    )

    return FilesAtVersionResponse(**response.json())


def delete_version_command(project_name: str, divbase_base_url: str, version_name: str) -> DeleteVersionResponse:
    """
    Delete a version from the bucket versioning file
    """
    request_data = DeleteVersionRequest(version_name=version_name)

    response = make_authenticated_request(
        method="DELETE",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/bucket-versions/delete?project_name={project_name}",
        json_data=request_data.model_dump(),
    )

    return DeleteVersionResponse(**response.json())


def list_files_command(divbase_base_url: str, project_name: str) -> list[str]:
    """List all files in a project."""
    response = make_authenticated_request(
        method="GET",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/s3/list?project_name={project_name}",
    )

    return response.json()


def download_files_command(
    divbase_base_url: str,
    project_name: str,
    all_files: list[str],
    download_dir: Path,
    bucket_version: str | None = None,
) -> list[Path]:
    """
    Download files from the given project's S3 bucket.
    """
    if not download_dir.is_dir():
        raise NotADirectoryError(
            f"The specified download directory '{download_dir}' is not a directory. Please create it or specify a valid directory before continuing."
        )

    # TODO - rewrite logic only once bucket versioning changes implemented for pre-signed url strategy
    if bucket_version:
        raise NotImplementedError("Downloading files at a specific bucket version is not yet (re)implemented.")

    query_params = {
        "project_name": project_name,
        "object_names": all_files,
    }

    response = make_authenticated_request(
        method="GET",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/s3/download?{urlencode(query_params, doseq=True)}",
    )

    return download_multiple_pre_signed_urls(pre_signed_urls=response.json(), download_dir=download_dir)


def upload_files_command(
    project_name: str, divbase_base_url: str, all_files: list[Path], safe_mode: bool
) -> dict[str, Path]:
    """
    Upload files to the project's S3 bucket.
    Files uploaded and there names in  returned as a list Paths

    Safe mode checks if any of the files that are to be uploaded already exist in the bucket.
    """
    object_names = [file.name for file in all_files]

    if safe_mode:
        response = make_authenticated_request(
            method="GET",
            divbase_base_url=divbase_base_url,
            api_route=f"v1/s3/list?project_name={project_name}",
        )
        current_files = response.json()

        existing_objects = set(object_names) & set(current_files)
        if existing_objects:
            raise FilesAlreadyInBucketError(existing_objects=list(existing_objects), project_name=project_name)

    query_params = {
        "project_name": project_name,
        "object_names": object_names,
    }

    response = make_authenticated_request(
        method="POST",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/s3/upload?{urlencode(query_params, doseq=True)}",
    )
    pre_signed_urls = response.json()

    return upload_multiple_pre_signed_urls(pre_signed_urls=pre_signed_urls, all_files=all_files)


def soft_delete_objects_command(divbase_base_url: str, project_name: str, all_files: list[str]) -> list[str]:
    """
    Soft delete objects from the project's S3 bucket.
    Returns a list of the soft deleted objects
    """
    query_params = {
        "project_name": project_name,
        "object_names": all_files,
    }

    response = make_authenticated_request(
        method="DELETE",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/s3/soft_delete?{urlencode(query_params, doseq=True)}",
    )
    return response.json().get("deleted", [])


def show_dimensions_command(project_config: ProjectConfig) -> dict[str, dict]:
    """
    Helper function used by the dimensions CLI command to show the dimensions index for a project.
    """
    s3_file_manager = create_s3_file_manager(project_config.s3_url)
    manager = VCFDimensionIndexManager(bucket_name=project_config.bucket_name, s3_file_manager=s3_file_manager)
    return manager.get_dimensions_info()
