"""
CLI commands for managing S3 bucket versions.
"""

from pathlib import Path

from divbase_tools.bucket_versioning import BucketVersionManager
from divbase_tools.exceptions import FilesAlreadyInBucketError, ObjectDoesNotExistInSpecifiedVersionError
from divbase_tools.s3_client import create_s3_file_manager
from divbase_tools.user_config import BucketConfig


def create_bucket_manager(bucket_config: BucketConfig) -> BucketVersionManager:
    """
    Helper function to create a BucketVersionManager instance.
    Used by the version and file subcommands of the CLI
    """
    s3_file_manager = create_s3_file_manager(bucket_config.s3_url)
    return BucketVersionManager(bucket_name=bucket_config.name, s3_file_manager=s3_file_manager)


def create_version_object_command(bucket_config: BucketConfig) -> None:
    manager = create_bucket_manager(bucket_config=bucket_config)
    manager.create_metadata_file()


def add_version_command(bucket_config: BucketConfig, name: str, description: str | None) -> None:
    manager = create_bucket_manager(bucket_config=bucket_config)
    manager.add_version(name=name, description=description)


def list_versions_command(bucket_config: BucketConfig) -> dict[str, dict]:
    manager = create_bucket_manager(bucket_config=bucket_config)
    return manager.get_version_info()


def list_files_at_version_command(bucket_config: BucketConfig, bucket_version: str) -> dict[str, str]:
    manager = create_bucket_manager(bucket_config=bucket_config)
    return manager.all_files_at_bucket_version(bucket_version=bucket_version)


def delete_version_command(bucket_config: BucketConfig, bucket_version: str) -> str:
    manager = create_bucket_manager(bucket_config=bucket_config)
    return manager.delete_version(bucket_version=bucket_version)


def list_files_command(bucket_config: BucketConfig) -> list[str]:
    s3_file_manager = create_s3_file_manager(url=bucket_config.s3_url)
    return s3_file_manager.list_files(bucket_name=bucket_config.name)


def download_files_command(
    bucket_config: BucketConfig, all_files: list[str], download_dir: Path, bucket_version: str | None = None
) -> list[Path]:
    s3_file_manager = create_s3_file_manager(url=bucket_config.s3_url)

    if not download_dir.is_dir():
        raise NotADirectoryError(
            f"The specified download directory '{download_dir}' is not a directory. Please create it or specify a valid directory before continuing."
        )

    # Get the files at the specified version.
    if bucket_version:
        bucket_version_manager = BucketVersionManager(bucket_name=bucket_config.name, s3_file_manager=s3_file_manager)
        files_at_version = bucket_version_manager.all_files_at_bucket_version(bucket_version=bucket_version)

        # check if all files specified exist for download exist at this bucket version
        missing_objects = [f for f in all_files if f not in files_at_version]
        if missing_objects:
            raise ObjectDoesNotExistInSpecifiedVersionError(
                bucket_name=bucket_config.name,
                bucket_version=bucket_version,
                missing_objects=missing_objects,
            )

        files_to_download = {file: files_at_version[file] for file in all_files}
    else:
        files_to_download = {file: None for file in all_files}

    download_files = s3_file_manager.download_files(
        objects=files_to_download,
        download_dir=download_dir,
        bucket_name=bucket_config.name,
    )

    return download_files


def upload_files_command(bucket_config: BucketConfig, all_files: list[Path], safe_mode: bool) -> dict[str, Path]:
    """
    Upload files to the specified S3 bucket.
    Files uploaded and there names in  returned as a list Paths

    Safe mode checks if any of the files that are to be uploaded already exist in the bucket.
    """
    s3_file_manager = create_s3_file_manager(url=bucket_config.s3_url)

    if safe_mode:
        bucket_version_manager = BucketVersionManager(bucket_name=bucket_config.name, s3_file_manager=s3_file_manager)
        current_files = bucket_version_manager._get_all_objects_names_and_ids().keys()
        file_names = [file.name for file in all_files]

        existing_objects = set(file_names) & set(current_files)
        if existing_objects:
            raise FilesAlreadyInBucketError(existing_objects=list(existing_objects), bucket_name=bucket_config.name)

    uploaded_files = s3_file_manager.upload_files(
        to_upload={file.name: file for file in all_files},
        bucket_name=bucket_config.name,
    )

    return uploaded_files


def delete_objects_command(bucket_config: BucketConfig, all_files: list[str]) -> list[str]:
    """
    Delete objects from the specified S3 bucket.
    Returns a list of the deleted objects
    """
    s3_file_manager = create_s3_file_manager(url=bucket_config.s3_url)
    return s3_file_manager.delete_objects(objects=all_files, bucket_name=bucket_config.name)
