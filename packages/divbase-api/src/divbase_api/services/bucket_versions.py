"""
Responsible for managing the versioning of the bucket state.

Note: this means the overall state of the bucket, not the individual files.

This was previously a part of the CLI but now due to redesign

"""

import logging
from dataclasses import dataclass
from datetime import datetime, timezone

import yaml

from divbase_api.config import settings

# TODO - move to API exceptions
from divbase_lib.exceptions import (
    BucketVersionAlreadyExistsError,
    BucketVersioningFileAlreadyExistsError,
    BucketVersionNotFoundError,
    ObjectDoesNotExistError,
)
from divbase_lib.s3_client import S3FileManager

VERSION_FILE_NAME = ".bucket_versions.yaml"

logger = logging.getLogger(__name__)


@dataclass
class BucketVersionManager:
    bucket_name: str
    s3_file_manager: S3FileManager
    version_info: dict

    def create_metadata_file(self, name: str | None, description: str | None) -> None:
        """
        Create the initial metadata file with a default version.
        """
        if self.version_info:
            raise BucketVersioningFileAlreadyExistsError(bucket_name=self.bucket_name)

        timestamp = self._create_timestamp()
        files = self._get_all_objects_names_and_ids()

        name = name or "v0.0.0"
        description = description or "First version"

        version_metadata = {"versions": {name: {"timestamp": timestamp, "description": description, "files": files}}}

        self._upload_bucket_version_file(version_data=version_metadata)
        logger.info("Bucket versioning file created and uploaded successfully.")

    def add_version(self, name: str, description: str | None) -> None:
        """
        Add a new version to the metadata file.
        """
        version_data = self.version_info
        if not version_data:  # aka file not created for this project yet.
            version_data["versions"] = {}

        if name in version_data["versions"]:
            raise BucketVersionAlreadyExistsError(bucket_name=self.bucket_name, version_name=name)

        timestamp = self._create_timestamp()
        files = self._get_all_objects_names_and_ids()
        _ = files.pop(VERSION_FILE_NAME, None)

        if not description:
            description = ""

        version_data["versions"][name] = {"timestamp": timestamp, "description": description, "files": files}
        self._upload_bucket_version_file(version_data=version_data)

    def delete_version(self, bucket_version: str) -> str:
        """
        Delete the version specified from the bucket versioning metadata file.

        Returns the version name deleted.
        """
        if not self.version_info:
            raise ObjectDoesNotExistError(key=VERSION_FILE_NAME, bucket_name=self.bucket_name)

        try:
            del self.version_info["versions"][bucket_version]
        except KeyError as err:
            raise BucketVersionNotFoundError(bucket_name=self.bucket_name, bucket_version=bucket_version) from err

        self._upload_bucket_version_file(version_data=self.version_info)
        return bucket_version

    def get_version_info(self) -> dict[str, dict]:
        if not self.version_info:
            logger.info("No bucket level versioning has been created for this bucket as of yet.")
            return {}
        return self.version_info["versions"]

    def all_files_at_bucket_version(self, bucket_version: str) -> dict[str, str]:
        """
        Get all files in the bucket at a user specified bucket version.
        If the version does not exist, raise an error.

        Return a dict of file names and their unique version IDs (hashes).
        """
        if not self.version_info:
            raise BucketVersionNotFoundError(bucket_name=self.bucket_name, bucket_version=bucket_version)

        try:
            version_info = self.version_info["versions"][bucket_version]
        except KeyError as err:
            raise BucketVersionNotFoundError(bucket_name=self.bucket_name, bucket_version=bucket_version) from err

        return version_info["files"]

    def _create_timestamp(self) -> str:
        return datetime.now(tz=timezone.utc).isoformat()

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
        self.s3_file_manager.upload_str_as_s3_object(
            key=VERSION_FILE_NAME, content=text_content, bucket_name=self.bucket_name
        )
        logging.info(f"New version updated in the bucket: {self.bucket_name}.")


def create_bucket_version_manager(bucket_name: str) -> BucketVersionManager:
    """
    Create a BucketVersionManager instance.

    If bucket version files exists, download version metadata file from the S3 bucket.
    If doesn't exist, create an empty version info dict.
    """
    s3_file_manager = S3FileManager(
        url=settings.s3.s3_internal_url,
        access_key=settings.s3.access_key.get_secret_value(),
        secret_key=settings.s3.secret_key.get_secret_value(),
    )

    try:
        content = s3_file_manager.download_s3_file_to_str(key=VERSION_FILE_NAME, bucket_name=bucket_name)
    except ObjectDoesNotExistError:
        logger.info(f"No bucket versioning file found in the bucket: {bucket_name}.")
        version_info = {}
    if not content:
        version_info = {}

    version_info = yaml.safe_load(content)
    return BucketVersionManager(bucket_name=bucket_name, s3_file_manager=s3_file_manager, version_info=version_info)
