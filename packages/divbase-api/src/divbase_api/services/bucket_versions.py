"""
Responsible for managing the versioning of the bucket state.

Note: A version entry here refers to the overall state of all files in the bucket at a timestamp.
So it is a record of all files and their unique version IDs (hashes) at a point in time.
"""

import logging
from dataclasses import dataclass
from datetime import datetime, timezone

import yaml

from divbase_api.api_config import settings
from divbase_api.exceptions import (
    BucketVersionAlreadyExistsError,
    BucketVersioningFileAlreadyExistsError,
    BucketVersionNotFoundError,
)
from divbase_api.services.s3_client import S3FileManager
from divbase_lib.api_schemas.bucket_versions import (
    AddVersionResponse,
    BucketVersionDetail,
    CreateVersioningFileResponse,
    FilesAtVersionResponse,
    VersionListResponse,
)
from divbase_lib.exceptions import ObjectDoesNotExistError

VERSION_FILE_NAME = ".bucket_versions.yaml"

logger = logging.getLogger(__name__)


@dataclass
class BucketVersionManager:
    bucket_name: str
    s3_file_manager: S3FileManager
    version_info: dict

    def create_metadata_file(self, name: str, description: str) -> CreateVersioningFileResponse:
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

        return CreateVersioningFileResponse(name=name, description=description)

    def add_version(self, name: str, description: str) -> AddVersionResponse:
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

        return AddVersionResponse(name=name, description=description)

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

    def get_version_info(self) -> VersionListResponse:
        if not self.version_info:
            logger.info("No bucket level versioning has been created for this bucket as of yet.")
            return VersionListResponse(versions={})

        # Convert raw version data to BucketVersionDetail objects
        versions = {}
        for name, data in self.version_info["versions"].items():
            versions[name] = BucketVersionDetail(
                timestamp=data["timestamp"], description=data["description"], files=data["files"]
            )

        return VersionListResponse(versions=versions)

    def all_files_at_bucket_version(self, bucket_version: str) -> FilesAtVersionResponse:
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

        return FilesAtVersionResponse(files=version_info["files"])

    def _create_timestamp(self) -> str:
        return datetime.now(tz=timezone.utc).isoformat()

    def _get_all_objects_names_and_ids(self) -> dict[str, str]:
        """
        Create a dict of the latest version of each file in the bucket and its unique versionID (hash).
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
        url=settings.s3.endpoint_url,
        access_key=settings.s3.access_key.get_secret_value(),
        secret_key=settings.s3.secret_key.get_secret_value(),
    )

    try:
        content = s3_file_manager.download_s3_file_to_str(key=VERSION_FILE_NAME, bucket_name=bucket_name)
        version_info = yaml.safe_load(content) if content else {}
    except ObjectDoesNotExistError:
        logger.info(f"No bucket versioning file found in the bucket: {bucket_name}.")
        version_info = {}

    return BucketVersionManager(bucket_name=bucket_name, s3_file_manager=s3_file_manager, version_info=version_info)
