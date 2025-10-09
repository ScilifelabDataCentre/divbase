"""
CLI commands for managing S3 bucket versions.
"""

from pathlib import Path
from urllib.parse import urlencode

from divbase_cli.bucket_versioning import BucketVersionManager
from divbase_cli.pre_signed_urls import download_multiple_pre_signed_urls, upload_multiple_pre_signed_urls
from divbase_cli.user_auth import make_authenticated_request
from divbase_cli.user_config import ProjectConfig
from divbase_lib.s3_client import create_s3_file_manager
from divbase_lib.vcf_dimension_indexing import VCFDimensionIndexManager


def create_bucket_manager(project_config: ProjectConfig) -> BucketVersionManager:
    """
    Helper function to create a BucketVersionManager instance.
    Used by the version and file subcommands of the CLI
    """
    s3_file_manager = create_s3_file_manager(project_config.s3_url)
    return BucketVersionManager(bucket_name=project_config.bucket_name, s3_file_manager=s3_file_manager)


def create_version_object_command(project_config: ProjectConfig) -> None:
    manager = create_bucket_manager(project_config=project_config)
    manager.create_metadata_file()


def add_version_command(project_config: ProjectConfig, name: str, description: str | None) -> None:
    manager = create_bucket_manager(project_config=project_config)
    manager.add_version(name=name, description=description)


def list_versions_command(project_config: ProjectConfig) -> dict[str, dict]:
    manager = create_bucket_manager(project_config=project_config)
    return manager.get_version_info()


def list_files_at_version_command(project_config: ProjectConfig, bucket_version: str) -> dict[str, str]:
    manager = create_bucket_manager(project_config=project_config)
    return manager.all_files_at_bucket_version(bucket_version=bucket_version)


def delete_version_command(project_config: ProjectConfig, bucket_version: str) -> str:
    manager = create_bucket_manager(project_config=project_config)
    return manager.delete_version(bucket_version=bucket_version)


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
    # TODO - reimplement safe mode after changes
    # Probably can do this by running list_files first?
    if safe_mode:
        raise NotImplementedError("Safe mode is not yet (re)implemented.")
    #     bucket_version_manager = BucketVersionManager(
    #         bucket_name=project_config.bucket_name, s3_file_manager=s3_file_manager
    #     )
    #     current_files = bucket_version_manager._get_all_objects_names_and_ids().keys()
    #     file_names = [file.name for file in all_files]

    #     existing_objects = set(file_names) & set(current_files)
    #     if existing_objects:
    #         raise FilesAlreadyInBucketError(
    #             existing_objects=list(existing_objects), bucket_name=project_config.bucket_name
    #         )

    object_names = [file.name for file in all_files]
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
