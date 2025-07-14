"""
An s3FileManager object that lets you do basic operations with a bucket.
"""

import os
from pathlib import Path

import boto3
import botocore

from divbase_tools.exceptions import DivBaseCredentialsNotFoundError, ObjectDoesNotExistError


class S3FileManager:
    def __init__(self, url: str, access_key: str, secret_key: str):
        self.s3_client = boto3.client(
            "s3",
            endpoint_url=url,
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

    def download_files(self, objects: dict[str, str | None], download_dir: Path, bucket_name: str) -> list[Path]:
        """
        Download objects/files from the S3 bucket to a local directory.

        objects is a dict with keys being the S3 object keys and values being the version IDs or None.
        """
        downloaded_files = []

        for key, version_id in objects.items():
            dest = download_dir / Path(key).name
            self._download_single_file(key, dest, bucket_name, version_id)
            downloaded_files.append(dest)

        return downloaded_files

    def download_s3_file_to_str(self, key: str, bucket_name: str) -> str:
        """
        Get the contents of a file from the S3 bucket as a string.
        """
        try:
            response = self.s3_client.get_object(Bucket=bucket_name, Key=key)
        except self.s3_client.exceptions.NoSuchKey as err:
            raise ObjectDoesNotExistError(key=key, bucket_name=bucket_name) from err
        return response["Body"].read().decode("utf-8")

    def upload_files(self, to_upload: dict[str, Path], bucket_name: str) -> dict[str, Path]:
        """
        Upload a list of files to the S3 bucket

        to_upload is a dict where keys are the name of the (to be) object and values is the path to the file to upload.
        Returns a dict with the keys and the source paths of the uploaded files.
        """
        uploaded_files = {}
        for key, source in to_upload.items():
            self.s3_client.upload_file(
                Filename=str(source),
                Bucket=bucket_name,
                Key=key,
            )
            uploaded_files[key] = source.resolve()

        return uploaded_files

    def upload_str_as_s3_object(self, key: str, content: str, bucket_name: str) -> None:
        """
        Upload a string (e.g. output of yaml.safe_dump()) as a new file to the S3 bucket.
        """
        self.s3_client.put_object(Bucket=bucket_name, Key=key, Body=content.encode("utf-8"))

    def delete_objects(self, objects: list[str], bucket_name: str) -> list[str]:
        """
        Delete files/objects from the S3 bucket, returns a list of the deleted objects.

        NOTE:
            - If the object was previously deleted (aka has a deletion marker) OR
            - If the object does not exist in the bucket,
        it will still be returned in the response as deleted and you cant distinguish between the two cases.
        Would need to check if the files exist to prevent this behaviour, but have to decide if this is worth doing.
        """
        delete_dict = {"Objects": [{"Key": obj} for obj in objects]}
        response = self.s3_client.delete_objects(
            Bucket=bucket_name,
            Delete=delete_dict,
        )
        deleted_object_names = [object_dict["Key"] for object_dict in response["Deleted"]]
        return deleted_object_names

    def delete_specific_object_versions(self, versioned_objects: dict[str, str], bucket_name: str) -> None:
        """Delete files from the S3 bucket by specifying the version of the object to delete."""
        delete_dict = {"Objects": []}
        for key, version in versioned_objects.items():
            delete_dict["Objects"].append({"Key": key, "VersionId": version})

        self.s3_client.delete_objects(
            Bucket=bucket_name,
            Delete=delete_dict,
        )

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

    def _download_single_file(self, key: str, dest: Path, bucket_name: str, version_id: str | None = None) -> str:
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


def create_s3_file_manager(
    url: str,
    s3_env_access_key_name: str = "DIVBASE_S3_ACCESS_KEY",
    s3_env_secret_key_name: str = "DIVBASE_S3_SECRET_KEY",
) -> S3FileManager:
    """
    Creates an S3FileManager instance using credentials from environment variables
    """
    access_key = os.getenv(s3_env_access_key_name)
    secret_key = os.getenv(s3_env_secret_key_name)

    if not access_key or not secret_key:
        raise DivBaseCredentialsNotFoundError(
            access_key_name=s3_env_access_key_name, secret_key_name=s3_env_secret_key_name
        )

    return S3FileManager(
        url=url,
        access_key=access_key,
        secret_key=secret_key,
    )
