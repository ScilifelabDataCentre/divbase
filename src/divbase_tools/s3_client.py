"""
An s3FileManager object that lets you do basic operations with a bucket.
"""

import os
from pathlib import Path

import boto3
import botocore
from dotenv import load_dotenv

from divbase_tools.exceptions import DivBaseCredentialsNotFoundError, ObjectDoesNotExistError
from divbase_tools.user_config import load_user_config

MINIO_URL = "api.divbase-testground.scilifelab-2-dev.sys.kth.se"


class S3FileManager:
    def __init__(self, url: str, access_key: str, secret_key: str):
        self.s3_client = boto3.client(
            "s3",
            endpoint_url=f"https://{url}",
            aws_access_key_id=access_key,
            aws_secret_access_key=secret_key,
        )

    def list_files(self, bucket_name: str) -> list[str]:
        """
        Return a list of all files in the S3 bucket.
        """
        files = []
        paginator = self.s3_client.get_paginator("list_objects_v2")
        for page in paginator.paginate(Bucket=bucket_name):
            for obj in page.get("Contents", []):
                files.append(obj["Key"])
        return files

    def download_file(self, key: str, dest: Path, bucket_name: str, version_id: str | None = None) -> str:
        """
        Download a file from S3 to a local path.
        Downloads the latest version of the file by default unless the the version_id is provided.

        Returns the key of the downloaded file.
        """
        extra_args = {"VersionId": version_id} if version_id else None
        try:
            self.s3_client.download_file(Bucket=bucket_name, Key=key, Filename=str(dest), ExtraArgs=extra_args)
        except botocore.exceptions.ClientError as err:
            if err.response["Error"]["Code"] == "404":
                raise ObjectDoesNotExistError(
                    key=key,
                    bucket_name=bucket_name,
                ) from err

            raise err

        return key

    def download_s3_file_to_str(self, key: str, bucket_name: str) -> str:
        """
        Get the contents of a file from the S3 bucket as a string.
        """
        try:
            response = self.s3_client.get_object(Bucket=bucket_name, Key=key)
        except self.s3_client.exceptions.NoSuchKey as err:
            raise ObjectDoesNotExistError(key=key, bucket_name=bucket_name) from err
        return response["Body"].read().decode("utf-8")

    def upload_file(self, key: str, source: Path, bucket_name: str) -> None:
        """
        Upload a file to the S3 bucket
        """
        self.s3_client.upload_file(
            Filename=str(source),
            Bucket=bucket_name,
            Key=key,
        )

    def upload_str_as_s3_object(self, key: str, content: str, bucket_name: str) -> None:
        """
        Upload a string (e.g. output of yaml.safe_dump()) as a new file to the S3 bucket.
        """
        self.s3_client.put_object(Bucket=bucket_name, Key=key, Body=content.encode("utf-8"))

    def latest_version_of_all_files(self, bucket_name: str) -> dict[str, str]:
        """
        Identify the latest version of each file in the bucket.
        Returns a dictionary with the file name as the key and the version ID as the value.
        """
        files = {}
        paginator = self.s3_client.get_paginator("list_object_versions")
        for page in paginator.paginate(Bucket=bucket_name):
            for obj in page.get("Versions", []):
                if obj["IsLatest"]:
                    files[obj["Key"]] = obj["VersionId"]
        return files


def get_credentials(access_key_name: str, secret_key_name: str) -> tuple[str, str]:
    """
    Get the S3 credentials from environment variables, unless they are provided as arguments.
    """
    load_dotenv()
    access_key = os.getenv(access_key_name)
    secret_key = os.getenv(secret_key_name)
    return access_key, secret_key


def config_to_s3_file_manager(config_path: Path) -> S3FileManager:
    """
    Creates an S3FileManager instance from the user config file.
    """
    user_config = load_user_config(config_path)

    access_key_name = user_config.get("DivBase_Access_Key_Env_Name")
    secret_key_name = user_config.get("DivBase_Secret_Key_Env_Name")

    access_key, secret_key = get_credentials(access_key_name=access_key_name, secret_key_name=secret_key_name)

    if not access_key or not secret_key:
        raise DivBaseCredentialsNotFoundError(
            access_key_name=access_key_name, secret_key_name=secret_key_name, config_path=config_path
        )

    return S3FileManager(
        url=MINIO_URL,
        access_key=access_key,
        secret_key=secret_key,
    )
