"""
CLI commands for managing S3 bucket versions.
"""

from pathlib import Path

from divbase_tools.bucket_versioning import BucketVersionManager
from divbase_tools.s3_client import S3FileManager


def create_version_object_command(bucket_name: str, s3_file_manager: S3FileManager) -> None:
    manager = BucketVersionManager(bucket_name=bucket_name, s3_file_manager=s3_file_manager)
    manager.create_metadata_file()


def add_version_command(bucket_name: str, name: str, description: str | None, s3_file_manager: S3FileManager) -> None:
    manager = BucketVersionManager(bucket_name=bucket_name, s3_file_manager=s3_file_manager)
    manager.add_version(name=name, description=description)


def list_versions_command(bucket_name: str, s3_file_manager: S3FileManager) -> None:
    manager = BucketVersionManager(bucket_name=bucket_name, s3_file_manager=s3_file_manager)
    manager.list_versions()


def delete_version_command(bucket_name: str, version_id: str, s3_file_manager: S3FileManager) -> None:
    # TODO should this even be a command?
    pass


def download_files_command(
    bucket_name: str, files: list[str], download_dir: Path, bucket_version: str, s3_file_manager: S3FileManager
) -> None:
    manager = BucketVersionManager(bucket_name=bucket_name, s3_file_manager=s3_file_manager)
    manager.download_files(files=files, download_dir=download_dir, bucket_version=bucket_version)


def upload_files_command(bucket_name: str, s3_file_manager: S3FileManager) -> None:
    pass
