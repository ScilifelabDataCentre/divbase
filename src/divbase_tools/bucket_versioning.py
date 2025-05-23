"""
Responsible for managing the versioning of the bucket state.

Note: this means the overall state of the bucket, not the individual files.
"""

from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import botocore
import yaml

from divbase_tools.s3_client import S3FileManager

VERSION_FILE_NAME = ".bucket_versions.yaml"


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
            print(f"Version file found inside Bucket: {self.bucket_name}.")

    def create_metadata_file(self):
        """Create the initial metadata file with a default version."""
        if self.version_info:
            print(f"Can't create a new version file as one already exists in the bucket: {self.bucket_name}.")
            return

        timestamp = self._create_timestamp()
        files = self._get_current_object_versions_ids()

        version_metadata = {
            "versions": {"v0.0.0": {"timestamp": timestamp, "description": "First version", "files": files}}
        }

        self._upload_bucket_version_file(version_data=version_metadata)
        print("Bucket versioning file created successfully.")

    def add_version(self, name: str, description: str | None) -> None:
        """Add a new version to the metadata file."""
        version_data = self.version_info

        timestamp = self._create_timestamp()
        files = self._get_current_object_versions_ids()

        if not description:
            description = ""
        version_data["versions"][name] = {"timestamp": timestamp, "description": description, "files": files}

        self._upload_bucket_version_file(version_data=version_data)
        # only update the version_info if the file was successfully written
        self.version_info = version_data

    def list_versions(self):
        if not self.version_info:
            print("No bucket level versioning has been created for this bucket as of yet.")
            return

        print("Bucket versions:")
        # TODO - perhaps user friendly timestamp format too?
        for version, details in self.version_info["versions"].items():
            if details["description"]:
                print(f"- {version}: {details['timestamp']} ({details['description']})")
            else:
                print(f"- {version}: {details['timestamp']} (No description provided)")

    def download_files(self, files: str, download_dir: str, bucket_version: str | None) -> None:
        """
        Given a list of comma separated files, download the files from the bucket into a specified download dir.
        """
        files_list = files.split(",")

        if bucket_version:
            try:
                version_info = self.version_info["versions"][bucket_version]
                print(f"Downloading the files: {files_list} from version: {bucket_version}")
            except KeyError:
                print(f"Version specified: {bucket_version} was not found in the bucket: {self.bucket_name}.")
                return
        else:
            print(f"Downloading the latest state of the files: {files_list}")
            version_info = None

        for file in files_list:
            download_path = Path(download_dir) / file

            if version_info:
                version_id = version_info["files"][file]
            else:
                version_id = None

            self.s3_file_manager.download_file(
                key=file, dest=download_path, bucket_name=self.bucket_name, version_id=version_id
            )

    def _create_timestamp(self) -> str:
        return datetime.now(tz=timezone.utc).isoformat()

    def _get_bucket_version_file(self) -> dict:
        """
        If bucket version files exists, download version metadata file from the S3 bucket.
        If doesn't exist, return empty dict.
        """
        content = self.s3_file_manager.download_s3_file_to_str(key=VERSION_FILE_NAME, bucket_name=self.bucket_name)
        if not content:
            return {}

        return yaml.safe_load(content)

    def _get_current_object_versions_ids(self) -> dict[str, str]:
        """
        Create a dict of the latest version of each file in the bucket and its unique versionID (hash).
        # TODO - does this handle soft deleted objects correctly?
        """
        files = self.s3_file_manager.latest_version_of_all_files()
        if not files:
            print(f"No files found in bucket '{self.bucket_name}'.")

        return files

    def _upload_bucket_version_file(self, version_data: dict) -> None:
        """
        Upload a new version of the bucket version metadata file to the S3 bucket.
        Works for both creating and updating the file.
        """
        text_content = yaml.safe_dump(version_data, sort_keys=False)
        try:
            self.s3_file_manager.upload_str_as_s3_object(key=VERSION_FILE_NAME, content=text_content)
            print(f"New version updated in the bucket: {self.bucket_name}.")
        except botocore.exceptions.ClientError as e:
            print(f"Failed to upload version file: {e}")
