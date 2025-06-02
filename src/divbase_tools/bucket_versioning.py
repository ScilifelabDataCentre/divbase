"""
Responsible for managing the versioning of the bucket state.

Note: this means the overall state of the bucket, not the individual files.
"""

import logging
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import botocore
import yaml

from divbase_tools.exceptions import (
    BucketVersionNotFoundError,
    FilesAlreadyInBucketError,
    ObjectDoesNotExistError,
    ObjectDoesNotExistInSpecifiedVersionError,
)
from divbase_tools.s3_client import S3FileManager

VERSION_FILE_NAME = ".bucket_versions.yaml"

logger = logging.getLogger(__name__)


@dataclass
class BucketVersionManager:
    bucket_name: str
    s3_file_manager: S3FileManager

    # not defined till post_init method run.
    version_info: dict = field(init=False)

    def __post_init__(self):
        # empty dict if file doesn't exist
        self.version_info = self._get_bucket_version_file()

        if self.version_info:
            logger.info(f"Version file found inside Bucket: {self.bucket_name}.")

    def create_metadata_file(self) -> None:
        """
        Create the initial metadata file with a default version.
        """
        if self.version_info:
            logger.error(f"Can't create a new version file as one already exists in the bucket: {self.bucket_name}.")
            return

        timestamp = self._create_timestamp()
        files = self._get_all_objects_names_and_ids()

        version_metadata = {
            "versions": {"v0.0.0": {"timestamp": timestamp, "description": "First version", "files": files}}
        }

        self._upload_bucket_version_file(version_data=version_metadata)
        logger.info("Bucket versioning file created and uploaded successfully.")

    def add_version(self, name: str, description: str | None) -> None:
        """
        Add a new version to the metadata file.
        """
        version_data = self.version_info

        timestamp = self._create_timestamp()
        files = self._get_all_objects_names_and_ids()
        _ = files.pop(VERSION_FILE_NAME, None)

        if not description:
            description = ""
        version_data["versions"][name] = {"timestamp": timestamp, "description": description, "files": files}

        self._upload_bucket_version_file(version_data=version_data)

    def get_version_info(self) -> dict[str, dict]:
        if not self.version_info:
            logger.info("No bucket level versioning has been created for this bucket as of yet.")
            return {}
        return self.version_info["versions"]

    def list_files_in_bucket(self) -> list[str]:
        file_list = self.s3_file_manager.list_files(bucket_name=self.bucket_name)

        if not file_list:
            logger.warning(f"No files found in bucket '{self.bucket_name}'.")

        return file_list

    def download_files(self, files: list[str], download_dir: str, bucket_version: str | None) -> list[str]:
        """
        Given a list of comma separated files, download the files from the bucket into a specified download dir.

        Return the files that were downloaded.
        """
        version_info = None
        if bucket_version:
            try:
                version_info = self.version_info["versions"][bucket_version]
                logger.info(f"Downloading files from bucket at version: {bucket_version}")
            except KeyError as err:
                logger.error(f"Version specified: {bucket_version} was not found in the bucket: {self.bucket_name}.")
                raise BucketVersionNotFoundError(bucket_name=self.bucket_name, bucket_version=bucket_version) from err

            # Validate all files specified exist for the bucket_version specified.
            missing_objects = [f for f in files if f not in version_info["files"]]
            if missing_objects:
                raise ObjectDoesNotExistInSpecifiedVersionError(
                    bucket_name=self.bucket_name,
                    bucket_version=bucket_version,
                    missing_objects=missing_objects,
                )

        downloaded_files = []
        for file in files:
            download_path = Path(download_dir) / file

            if version_info:
                version_id = version_info["files"][file]
            else:
                version_id = None

            dloaded_file = self.s3_file_manager.download_file(
                key=file, dest=download_path, bucket_name=self.bucket_name, version_id=version_id
            )
            downloaded_files.append(dloaded_file)
        return downloaded_files

    def upload_files(self, files: list[Path]) -> list[Path]:
        """
        Upload files to the bucket.
        """
        uploaded_files = []
        for file in files:
            file_name = file.name
            self.s3_file_manager.upload_file(key=file_name, source=file, bucket_name=self.bucket_name)
            uploaded_files.append(file)
            logger.info(f"Uploaded file: {file.resolve()} to bucket: {self.bucket_name} as {file_name}.")
        return uploaded_files

    def safe_upload_files(self, files: list[Path]) -> list[Path]:
        """
        Upload files to the bucket.
        First check if any of the files already exist in the bucket, exit early if do.
        """
        current_files = self._get_all_objects_names_and_ids().keys()
        file_names = [file.name for file in files]

        existing_objects = set(file_names) & set(current_files)
        if existing_objects:
            raise FilesAlreadyInBucketError(existing_objects=list(existing_objects), bucket_name=self.bucket_name)

        return self.upload_files(files=files)

    def _create_timestamp(self) -> str:
        return datetime.now(tz=timezone.utc).isoformat()

    def _get_bucket_version_file(self) -> dict:
        """
        If bucket version files exists, download version metadata file from the S3 bucket.
        If doesn't exist, return empty dict.
        """
        try:
            content = self.s3_file_manager.download_s3_file_to_str(key=VERSION_FILE_NAME, bucket_name=self.bucket_name)
        except ObjectDoesNotExistError:
            logger.info(f"No bucket versioning file found in the bucket: {self.bucket_name}.")
            return {}
        if not content:
            return {}

        return yaml.safe_load(content)

    def _get_all_objects_names_and_ids(self) -> dict[str, str]:
        """
        Create a dict of the latest version of each file in the bucket and its unique versionID (hash).
        # TODO - does this handle soft deleted objects correctly?
        """
        files = self.s3_file_manager.latest_version_of_all_files(bucket_name=self.bucket_name)
        if not files:
            logging.info(f"No files found in bucket '{self.bucket_name}'.")
        return files

    def _upload_bucket_version_file(self, version_data: dict) -> None:
        """
        Upload a new version of the bucket version metadata file to the S3 bucket.
        Works for both creating and updating the file.
        """
        text_content = yaml.safe_dump(version_data, sort_keys=False)
        try:
            self.s3_file_manager.upload_str_as_s3_object(
                key=VERSION_FILE_NAME, content=text_content, bucket_name=self.bucket_name
            )
            logging.info(f"New version updated in the bucket: {self.bucket_name}.")
        except botocore.exceptions.ClientError as e:
            logging.error(f"Failed to upload bucket version file: {e}")
