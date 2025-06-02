"""
CLI commands for managing S3 bucket versions.
"""

from pathlib import Path

from divbase_tools.bucket_versioning import BucketVersionManager
from divbase_tools.s3_client import config_to_s3_file_manager


def create_bucket_manager(bucket_name: str, config_path: Path) -> BucketVersionManager:
    """
    Helper function to create a BucketVersionManager instance.
    Used by the version and file subcommands of the CLI
    """
    s3_file_manager = config_to_s3_file_manager(config_path=config_path)
    return BucketVersionManager(bucket_name=bucket_name, s3_file_manager=s3_file_manager)


def create_version_object_command(bucket_name: str, config_path: Path) -> None:
    manager = create_bucket_manager(bucket_name=bucket_name, config_path=config_path)
    manager.create_metadata_file()


def add_version_command(bucket_name: str, name: str, description: str | None, config_path: Path) -> None:
    manager = create_bucket_manager(bucket_name=bucket_name, config_path=config_path)
    manager.add_version(name=name, description=description)


def list_versions_command(bucket_name: str, config_path: Path) -> dict[str, dict]:
    manager = create_bucket_manager(bucket_name=bucket_name, config_path=config_path)
    return manager.get_version_info()


def delete_version_command(bucket_name: str, version_name: str, config_path: Path) -> None:
    # TODO should this even be a command?
    pass


def list_files_command(bucket_name: str, config_path: Path) -> None:
    manager = create_bucket_manager(bucket_name=bucket_name, config_path=config_path)
    return manager.list_files_in_bucket()


def download_files_command(
    bucket_name: str, all_files: list[str], download_dir: Path, bucket_version: str, config_path: Path
) -> list[str]:
    manager = create_bucket_manager(bucket_name=bucket_name, config_path=config_path)
    downloaded_files = manager.download_files(files=all_files, download_dir=download_dir, bucket_version=bucket_version)
    return downloaded_files


def upload_files_command(bucket_name: str, all_files: list[Path], safe_mode: bool, config_path: Path) -> list[Path]:
    """
    Upload files to the specified S3 bucket.
    Files uploaded returned as a list Paths
    """
    manager = create_bucket_manager(bucket_name=bucket_name, config_path=config_path)

    if safe_mode:
        uploaded_files = manager.safe_upload_files(files=all_files)
    else:
        uploaded_files = manager.upload_files(files=all_files)

    return uploaded_files


def delete_files_command(bucket_name: str, config_path: Path) -> None:
    # TODO should this even be a command?
    pass
