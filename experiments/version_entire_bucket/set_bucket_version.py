"""
Manage versions of the entire bucket

This script is interacted with via argparse.

Versioning specified by updating a file in the bucket called ".bucket_versions.yaml".

Versioning done by setting the version according to the current timestamp.
The version file is always kept in the bucket, not stored on disk.
"""

import argparse
import os
from dataclasses import dataclass, field
from datetime import datetime, timezone

import boto3
import botocore
import yaml
from botocore.exceptions import ClientError
from dotenv import load_dotenv

VERSION_FILE_NAME = ".bucket_versions.yaml"
MINIO_URL = "api.divbase-testground.scilifelab-2-dev.sys.kth.se"

ACCESS_KEY = os.getenv("MINIO_SQUIRREL_USER_ACCESS_KEY")
SECRET_KEY = os.getenv("MINIO_SQUIRREL_USER_SECRET_KEY")


@dataclass
class BucketVersionManager:
    bucket_name: str
    # not defined till post_init method run.
    s3_client: "botocore.client.S3" = field(init=False)
    version_info: dict = field(init=False)

    def __post_init__(self):
        load_dotenv()
        self.s3_client = boto3.client(
            "s3",
            endpoint_url=f"https://{MINIO_URL}",
            aws_access_key_id=ACCESS_KEY,
            aws_secret_access_key=SECRET_KEY,
        )

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
        version_data = {"versions": {"v0.1.0": {"timestamp": timestamp, "description": "Initial version"}}}

        self._upload_bucket_version_file(version_data=version_data)
        print("Bucket versioning file created successfully.")

    def add_version(self, name: str, description: str | None) -> None:
        """Add a new version to the metadata file."""
        version_data = self.version_info

        timestamp = self._create_timestamp()
        if not description:
            description = ""
        version_data["versions"][name] = {"timestamp": timestamp, "description": description}

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
            print(f"- {version}: {details['timestamp']} ({details['description']})")

    def _create_timestamp(self) -> str:
        return datetime.now(tz=timezone.utc).isoformat()

    def _get_bucket_version_file(self) -> dict:
        """
        If bucket version files exists, download version metadata file exists in the S3 bucket.
        If doesn't exist, return empty dict.
        """
        try:
            response = self.s3_client.get_object(Bucket=self.bucket_name, Key=VERSION_FILE_NAME)
        except self.s3_client.exceptions.NoSuchKey:
            print(f"The version file was not found in bucket: {self.bucket_name}, hopefully this is correct...")
            return {}

        content = response["Body"].read().decode("utf-8")
        # TODO - how would this handle an empty file, aka file exists but is empty?
        return yaml.safe_load(content)

    def _upload_bucket_version_file(self, version_data: dict) -> None:
        """
        Upload a version metadata file to the S3 bucket.
        Works for both creating and updating the file.
        """
        try:
            self.s3_client.put_object(Bucket=self.bucket_name, Key=VERSION_FILE_NAME, Body=yaml.safe_dump(version_data))
            print(f"New version updated in the bucket: {self.bucket_name}.")
        except ClientError as e:
            print(f"Failed to upload version file: {e}")


def create_command(args):
    manager = BucketVersionManager(args.bucket_name)
    manager.create_metadata_file()


def add_command(args):
    manager = BucketVersionManager(args.bucket_name)
    manager.add_version(args.name, args.description)


def list_command(args):
    manager = BucketVersionManager(args.bucket_name)
    manager.list_versions()


def delete_command(args):
    # TODO should this even be a command?
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Manage bucket versions")

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("--bucket_name", default=".", help="Name of the storage bucket for the project.")

    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: create
    create_parser = subparsers.add_parser(
        "create", parents=[parent_parser], help="Create the versioning metadata file."
    )
    create_parser.set_defaults(func=create_command)

    # Subcommand: add
    add_parser = subparsers.add_parser("add", parents=[parent_parser], help="Add a new version.")
    add_parser.add_argument("--name", required=True, help="Name of the version (e.g., semantic version).")
    add_parser.add_argument("--description", required=False, help="Optional description of the version.")
    add_parser.set_defaults(func=add_command)

    # Subcommand: list
    list_parser = subparsers.add_parser("list", parents=[parent_parser], help="List all versions.")
    list_parser.set_defaults(func=list_command)

    args = parser.parse_args()
    args.func(args)
