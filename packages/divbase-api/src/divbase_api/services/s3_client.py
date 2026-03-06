"""
An s3FileManager object that lets you do basic operations with a bucket.

This class is a wrapper around the boto3 S3 client.

See the 'docs/development/s3_transfers.md' for more info on how S3 transfers are handled.
"""

import logging
import os
from datetime import datetime
from pathlib import Path

import boto3
import stamina
from boto3.s3.transfer import TransferConfig
from botocore.config import Config

from divbase_api.exceptions import DownloadedFileChecksumMismatchError, ObjectDoesNotExistError
from divbase_lib.api_schemas.s3 import (
    ListObjectsResponse,
    ObjectDetails,
    ObjectInfoResponse,
    ObjectVersionInfo,
    RestoreObjectsResponse,
    SoftDeletedObjectDetails,
)
from divbase_lib.divbase_constants import S3_MULTIPART_CHUNK_SIZE, S3_MULTIPART_UPLOAD_THRESHOLD
from divbase_lib.exceptions import ChecksumVerificationError
from divbase_lib.s3_checksums import verify_downloaded_checksum

logger = logging.getLogger(__name__)

S3_BATCH_SIZE = 1000


def retry_on_retriable_checksum_errors(exception: Exception) -> bool:
    """
    Retry condition function for stamina to decorators to only retry on checksum verification errors.

    (We don't need retry logic for other S3 errors as boto3 has inbuilt retry logic for these.)
    """
    return isinstance(exception, DownloadedFileChecksumMismatchError)


class S3FileManager:
    """An S3 client wrapper to do basic S3 operations."""

    def __init__(self, url: str, access_key: str, secret_key: str):
        self.s3_client = boto3.client(
            "s3",
            endpoint_url=url,
            aws_access_key_id=access_key,
            aws_secret_access_key=secret_key,
            config=Config(
                retries={
                    "max_attempts": 5,
                    "mode": "adaptive",
                }
            ),
        )
        # multipart up/download config
        self.transfer_config = TransferConfig(
            multipart_threshold=S3_MULTIPART_UPLOAD_THRESHOLD,
            multipart_chunksize=S3_MULTIPART_CHUNK_SIZE,
            max_concurrency=10,
            use_threads=True,
        )

    def list_files(self, bucket_name: str) -> list[str]:
        """
        Return a list of all files in the S3 bucket.

        This will run multiple requests if there are more than S3_BATCH_SIZE files in the bucket.
        Used by worker to list all vcf.gz files for processing.
        """
        files = []
        paginator = self.s3_client.get_paginator("list_objects_v2")
        for page in paginator.paginate(Bucket=bucket_name):
            for obj in page.get("Contents", []):
                files.append(obj["Key"])
        return files

    def list_files_detailed(
        self, bucket_name: str, prefix: str | None = None, next_token: str | None = None
    ) -> ListObjectsResponse:
        """
        Return a list of up to S3_BATCH_SIZE files in the S3 bucket with detailed info about each file.
        This is used by CLI users via the API.

        Pagination is supported via the next_token parameter, so a client may need to make multiple calls to get all files.
        """
        request_args = {
            "Bucket": bucket_name,
            "MaxKeys": S3_BATCH_SIZE,
        }
        if prefix:
            request_args["Prefix"] = prefix
        if next_token:
            request_args["ContinuationToken"] = next_token

        response = self.s3_client.list_objects_v2(**request_args)

        items = []
        for obj in response.get("Contents", []):
            items.append(
                ObjectDetails(
                    name=obj["Key"],
                    size=obj["Size"],
                    last_modified=obj["LastModified"],
                    etag=obj["ETag"].strip('"'),
                )
            )

        new_next_token: str | None = response.get("NextContinuationToken")
        return ListObjectsResponse(objects=items, next_token=new_next_token)

    def list_soft_deleted_files(self, bucket_name: str) -> list[SoftDeletedObjectDetails]:
        """
        list all soft-deleted filesobjects in a bucket.
        A soft-deleted object is one whose latest S3 object version is a delete marker.

        NOTE: As we expect the number of soft-deleted files to be low, we can just paginate here,
        rather than have the client make multiple calls with next tokens like in 'list_files_detailed'.
        """
        paginator = self.s3_client.get_paginator("list_object_versions")
        soft_deleted_files = []

        for page in paginator.paginate(Bucket=bucket_name):
            for marker in page.get("DeleteMarkers", []):
                if marker.get("IsLatest"):
                    soft_deleted_files.append(
                        SoftDeletedObjectDetails(
                            name=marker["Key"],
                            last_modified=marker["LastModified"],
                        )
                    )
        return soft_deleted_files

    def get_detailed_object_info(self, bucket_name: str, object_name: str) -> ObjectInfoResponse:
        """
        Retrieves all versions of a specific object in an S3 bucket, with details about each version in the bucket.
        Included in the response object is whether the object is currently deleted (i.e., has a deletion marker as the latest version).

        This is the level of detail a user needs to know about (i.e., then don't need to know about older deletion markers).
        """
        object_exists = False
        is_currently_deleted = False

        response = self.s3_client.list_object_versions(Bucket=bucket_name, Prefix=object_name)

        all_versions = []
        for version in response.get("Versions", []):
            if version["Key"] == object_name:
                object_exists = True
                all_versions.append(
                    ObjectVersionInfo(
                        version_id=version["VersionId"],
                        is_latest=version["IsLatest"],
                        last_modified=version["LastModified"],
                        size=version["Size"],
                        etag=version["ETag"].strip('"'),
                    )
                )

        if not object_exists:
            raise ObjectDoesNotExistError(key=object_name, bucket_name=bucket_name)

        for delete_marker in response.get("DeleteMarkers", []):
            if delete_marker["Key"] == object_name and delete_marker["IsLatest"]:
                is_currently_deleted = True
                break

        # Sort the results object by last_modified date
        all_versions.sort(key=lambda v: v.last_modified, reverse=True)

        return ObjectInfoResponse(
            object_name=object_name,
            is_currently_deleted=is_currently_deleted,
            versions=all_versions,
        )

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
                Config=self.transfer_config,
            )
            uploaded_files[key] = source.resolve()

        return uploaded_files

    def soft_delete_objects(self, objects: list[str], bucket_name: str) -> list[str]:
        """
        Soft delete files/objects from the S3 bucket, returns a list of the deleted objects.

        NOTE:
            - If the object was previously deleted (aka has a deletion marker) OR
            - If the object does not exist in the bucket,
        it will still be returned in the response as deleted and you cant distinguish between the two cases.
        Would need to check if the files exist to prevent this behaviour, but have to decide if this is worth doing.
        """
        deleted_object_names = []
        for i in range(0, len(objects), S3_BATCH_SIZE):
            batch = objects[i : i + S3_BATCH_SIZE]
            delete_dict = {"Objects": [{"Key": obj} for obj in batch]}
            response = self.s3_client.delete_objects(
                Bucket=bucket_name,
                Delete=delete_dict,
            )
            deleted_object_names.extend([object_dict["Key"] for object_dict in response.get("Deleted", [])])
        return deleted_object_names

    def restore_objects(self, objects: list[str], bucket_name: str) -> RestoreObjectsResponse:
        """
        Attempts to restore a list of soft-deleted objects by removing their delete markers.

        This method categorizes the outcome for each object into one of two states:
        - 'restored':
            The object had a delete marker which was successfully removed. OR
            The object was already live (no delete marker present as latest version).
        - 'not_restored':
            The object does not exist in the bucket (no versions or markers) OR
            we failed to remove the delete marker due to an unexpected s3 error.
        """
        restored_objects, not_restored = [], []

        markers_to_delete = []
        for obj_key in objects:
            response = self.s3_client.list_object_versions(Bucket=bucket_name, Prefix=obj_key)

            # Filters to make sure working with the exact object key, not just objects with the same prefix
            delete_markers = [marker for marker in response.get("DeleteMarkers", []) if marker["Key"] == obj_key]
            versions = [version for version in response.get("Versions", []) if version["Key"] == obj_key]

            is_deleted = delete_markers and delete_markers[0]["IsLatest"]
            is_not_deleted = versions and versions[0]["IsLatest"]

            if is_deleted:
                version_id = delete_markers[0]["VersionId"]
                markers_to_delete.append({"Key": obj_key, "VersionId": version_id})
            elif is_not_deleted:
                restored_objects.append(obj_key)
            else:
                # object not found, could be typo or was already hard deleted.
                not_restored.append(obj_key)

        if not markers_to_delete:
            return RestoreObjectsResponse(restored=restored_objects, not_restored=not_restored)

        for i in range(0, len(markers_to_delete), S3_BATCH_SIZE):
            batch = markers_to_delete[i : i + S3_BATCH_SIZE]
            delete_response = self.s3_client.delete_objects(
                Bucket=bucket_name,
                Delete={"Objects": batch},
            )
            for deleted_obj in delete_response.get("Deleted", []):
                restored_objects.append(deleted_obj["Key"])

            if delete_response.get("Errors"):
                for error in delete_response["Errors"]:
                    logger.error(
                        f"Failed to remove delete marker for '{error['Key']}' (version: {error['VersionId']}): "
                        f"{error['Code']} - {error['Message']}"
                    )
                    not_restored.append(error["Key"])

        return RestoreObjectsResponse(restored=restored_objects, not_restored=not_restored)

    def hard_delete_specific_object_versions(self, objects: list[dict[str, str]], bucket_name: str) -> None:
        """
        Hard delete a file from the S3 bucket by specifying the version of the object to delete.

        The list of objects to delete should have the form:
        [
            {"Key": "file1.txt", "VersionId": "version-id-1"},
            {"Key": "file2.txt", "VersionId": "version-id-2"},
        ]

        For versioned buckets,
        - deleting an object without specifying a version id just adds a deletion marker.
        - If you specify the version id, it hard deletes that specific version of the object.
        """
        for i in range(0, len(objects), S3_BATCH_SIZE):
            batch = objects[i : i + S3_BATCH_SIZE]
            self.s3_client.delete_objects(Bucket=bucket_name, Delete={"Objects": batch})

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

    @stamina.retry(on=retry_on_retriable_checksum_errors, attempts=3)
    def _download_single_file(self, key: str, dest: Path, bucket_name: str, version_id: str | None = None) -> str:
        """
        Download a file from S3 to a local path.
        Downloads the latest version of the file by default unless the the version_id is provided.

        Returns the key of the downloaded file.

        A failed checksum verification will raise a DownloadedFileChecksumMismatchError, which will trigger a retry with stamina.
        boto3 has its own retry logic for other S3 errors.
        """
        extra_args = {"VersionId": version_id} if version_id else None
        try:
            self.s3_client.download_file(
                Bucket=bucket_name,
                Key=key,
                Filename=str(dest),
                ExtraArgs=extra_args,
                Config=self.transfer_config,
            )
        except self.s3_client.exceptions.ClientError as err:
            if err.response["Error"]["Code"] == "404":
                raise ObjectDoesNotExistError(
                    key=key,
                    bucket_name=bucket_name,
                ) from err

            raise err

        # validate the etag on s3 matches the downloaded file
        expected_checksum = self.get_object_checksum_if_exists(
            bucket_name=bucket_name,
            object_name=key,
            version_id=version_id,
        )
        if not expected_checksum:
            raise ObjectDoesNotExistError(key=key, bucket_name=bucket_name)

        try:
            verify_downloaded_checksum(file_path=dest, expected_checksum=expected_checksum)
        except ChecksumVerificationError as err:
            dest.unlink()  # delete the invalid file as will try again instead.
            logger.warning(
                f"Checksum verification failed for downloaded file '{dest}'. "
                f"Expected: {err.expected_checksum}, Calculated: {err.calculated_checksum}.  Retrying download..."
            )
            raise DownloadedFileChecksumMismatchError(
                file_path=dest,
                expected_checksum=err.expected_checksum,
                calculated_checksum=err.calculated_checksum,
            ) from None

        return key

    def get_object_checksum_if_exists(
        self, bucket_name: str, object_name: str, version_id: str | None = None
    ) -> str | None:
        """
        Get the MD5 checksum of a file in the bucket, if it exists.
        Returns None if the file does not exist.
        """
        args = {"Bucket": bucket_name, "Key": object_name}
        if version_id:
            args["VersionId"] = version_id
        try:
            response = self.s3_client.head_object(**args)
        except self.s3_client.exceptions.ClientError as err:
            if err.response["Error"]["Code"] == "404":
                return None
            else:
                logger.error(
                    f"Unexpected error getting checksum for object '{object_name}' in bucket '{bucket_name}': {err}"
                )
                return None
        return response.get("ETag", "").strip('"')

    def get_expired_soft_deleted_objects(
        self, bucket_name: str, cutoff_date: datetime, prefix_exclude: str | None = None
    ) -> set[str]:
        """
        Identify "expired" soft-deleted objects in a bucket.

        A soft-deleted object is one whose latest S3 object version is a delete marker.
        A soft-deleted object is considered "expired" if its delete marker's LastModified timestamp is earlier than the provided cutoff_date.
        """
        paginator = self.s3_client.get_paginator("list_object_versions")
        expired_objects = set()

        for page in paginator.paginate(Bucket=bucket_name):
            # Check for delete markers that are the "Latest" version
            for marker in page.get("DeleteMarkers", []):
                if marker.get("IsLatest") and marker["LastModified"] < cutoff_date:
                    object_key = marker["Key"]
                    if prefix_exclude and object_key.startswith(prefix_exclude):
                        continue
                    expired_objects.add(object_key)
        return expired_objects

    def get_bucket_usage_bytes(self, bucket_name: str) -> int:
        """
        Get the total size in bytes of all objects in the bucket.

        Unfortunately there is no single API call to get the total bucket size, hence the summing approach below.
        """
        total_size = 0
        paginator = self.s3_client.get_paginator("list_objects_v2")
        for page in paginator.paginate(Bucket=bucket_name):
            for obj in page.get("Contents", []):
                total_size += obj["Size"]
        return total_size


def create_s3_file_manager(url: str) -> S3FileManager:
    """Helper function to creates an S3FileManager instance using the S3 service account's credentials"""
    access_key = os.environ["S3_SERVICE_ACCOUNT_ACCESS_KEY"]
    secret_key = os.environ["S3_SERVICE_ACCOUNT_SECRET_KEY"]

    return S3FileManager(
        url=url,
        access_key=access_key,
        secret_key=secret_key,
    )
