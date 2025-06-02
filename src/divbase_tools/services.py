"""
CLI commands for managing S3 bucket versions.
"""

from pathlib import Path

from divbase_tools.bucket_versioning import BucketVersionManager
from divbase_tools.exceptions import FilesAlreadyInBucketError, ObjectDoesNotExistInSpecifiedVersionError
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
    s3_file_manager = config_to_s3_file_manager(config_path=config_path)

    # Get the files at the specified version.
    if bucket_version:
        bucket_version_manager = BucketVersionManager(bucket_name=bucket_name, s3_file_manager=s3_file_manager)
        files_at_version = bucket_version_manager.all_files_at_bucket_version(bucket_version=bucket_version)

        # check if all files specified exist for download exist a this bucket version
        missing_objects = [f for f in all_files if f not in files_at_version]
        if missing_objects:
            raise ObjectDoesNotExistInSpecifiedVersionError(
                bucket_name=bucket_name,
                bucket_version=bucket_version,
                missing_objects=missing_objects,
            )

        files_to_download = {file: files_at_version[file] for file in all_files}
    else:
        files_to_download = {file: None for file in all_files}

    download_files = []
    for file_name, version_id in files_to_download.items():
        destination_path = download_dir / file_name
        dloaded_file = s3_file_manager.download_file(
            key=file_name, dest=destination_path, bucket_name=bucket_name, version_id=version_id
        )
        download_files.append(dloaded_file)

    return download_files


def upload_files_command(bucket_name: str, all_files: list[Path], safe_mode: bool, config_path: Path) -> list[Path]:
    """
    Upload files to the specified S3 bucket.
    Files uploaded returned as a list Paths

    Safe mode checks if any of the files that are to be uploaded already exist in the bucket.
    """
    s3_file_manager = config_to_s3_file_manager(config_path=config_path)

    if safe_mode:
        bucket_version_manager = BucketVersionManager(bucket_name=bucket_name, s3_file_manager=s3_file_manager)
        current_files = bucket_version_manager._get_all_objects_names_and_ids().keys()
        file_names = [file.name for file in all_files]

        existing_objects = set(file_names) & set(current_files)
        if existing_objects:
            raise FilesAlreadyInBucketError(existing_objects=list(existing_objects), bucket_name=bucket_name)

    uploaded_files = []
    for file_path in all_files:
        file_name = file_path.name
        s3_file_manager.upload_file(key=file_name, source=file_path, bucket_name=bucket_name)
        uploaded_files.append(file_path.resolve())
    return uploaded_files


def delete_files_command(bucket_name: str, config_path: Path) -> None:
    # TODO should this even be a command?
    pass
