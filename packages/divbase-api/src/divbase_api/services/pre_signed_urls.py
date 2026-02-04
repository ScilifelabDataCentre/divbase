"""
Logic to generate presigned URLs for S3 objects.
Gives end users access to S3 objects without needing AWS credentials.

Certain operations like "list files" etc... instead performed directly by the backend and returned to end user.
Pre-signed url approach not used when the request to s3 can be done in the (user to API) request-response cycle.
"""

import logging
from functools import lru_cache
from math import ceil

import boto3
from botocore.config import Config

from divbase_api.api_config import settings
from divbase_lib.api_schemas.s3 import (
    AbortMultipartUploadResponse,
    CompleteMultipartUploadResponse,
    CreateMultipartUploadResponse,
    PreSignedDownloadResponse,
    PreSignedSinglePartUploadResponse,
    PresignedUploadPartUrlResponse,
    UploadedPart,
)
from divbase_lib.divbase_constants import (
    DOWNLOAD_URL_EXPIRATION_SECONDS,
    MULTI_PART_UPLOAD_URL_EXPIRATION_SECONDS,
    S3_MULTIPART_CHUNK_SIZE,
    SINGLE_PART_UPLOAD_URL_EXPIRATION_SECONDS,
)

logger = logging.getLogger(__name__)


class S3PreSignedService:
    """
    Service to create pre-signed URLs for S3 object upload/download.
    Knows nothing about user authentication/authorization.

    We have 2 s3 client instances here:
    one for presigning urls,
    one for direct s3 operations with the service account.
    Both in production and in local dev these point to different endpoint URLs.
    """

    def __init__(self):
        s3_config = Config(
            retries={
                "max_attempts": 5,
                "mode": "adaptive",
            }
        )
        self.s3_pre_signing_client = boto3.client(
            "s3",
            endpoint_url=settings.s3.presigning_url,
            aws_access_key_id=settings.s3.access_key.get_secret_value(),
            aws_secret_access_key=settings.s3.secret_key.get_secret_value(),
            config=s3_config,
        )
        self.s3_client = boto3.client(
            "s3",
            endpoint_url=settings.s3.endpoint_url,
            aws_access_key_id=settings.s3.access_key.get_secret_value(),
            aws_secret_access_key=settings.s3.secret_key.get_secret_value(),
            config=s3_config,
        )

    def create_presigned_url_for_download(
        self, bucket_name: str, object_name: str, version_id: str | None
    ) -> PreSignedDownloadResponse:
        """
        Generate a presigned URL for S3 object download.
        Client determines whether to do single-part or multi-part download from this single URL.

        The generate_presigned_url method from boto3 has following params:
            :param client_method_name: Name of the S3.Client method, e.g., 'list_buckets'
            :param method_parameters: Dictionary of parameters to send to the method
            :param expiration: Time in seconds for the presigned URL to remain valid
            :return: Presigned URL as string. If error, returns None.
        """
        extra_args = {"VersionId": version_id} if version_id else None
        url = self.s3_pre_signing_client.generate_presigned_url(
            ClientMethod="get_object",
            Params={"Bucket": bucket_name, "Key": object_name, **(extra_args or {})},
            ExpiresIn=DOWNLOAD_URL_EXPIRATION_SECONDS,
        )

        return PreSignedDownloadResponse(
            name=object_name,
            pre_signed_url=url,
            version_id=version_id,
        )

    def create_presigned_url_for_single_part_upload(
        self, bucket_name: str, object_name: str, content_length: int, md5_hash: str | None = None
    ) -> PreSignedSinglePartUploadResponse:
        """
        Generate a presigned S3 PUT URL to upload a file to S3 (single part upload).
        The response object contains the object name, pre-signed URL, and any headers that must be included in the PUT request.

        NOTE:
        - If the headers are not included then S3 will return a 403 Forbidden error due to signature mismatch.
        - An initial attempt to use pre-signed POST URLs was given up on due to NetApp not seeming to support them.
        """
        put_headers: dict[str, str] = {}
        upload_args: dict[str, str | int] = {"Bucket": bucket_name, "Key": object_name}

        upload_args["ContentType"] = "application/octet-stream"
        put_headers["Content-Type"] = "application/octet-stream"

        upload_args["ContentLength"] = content_length  # boto3 expects this as an int.
        put_headers["Content-Length"] = str(content_length)

        if md5_hash:
            upload_args["ContentMD5"] = md5_hash
            put_headers["Content-MD5"] = md5_hash

        pre_signed_url = self.s3_pre_signing_client.generate_presigned_url(
            HttpMethod="PUT",
            ClientMethod="put_object",
            Params=upload_args,
            ExpiresIn=SINGLE_PART_UPLOAD_URL_EXPIRATION_SECONDS,
        )
        return PreSignedSinglePartUploadResponse(
            name=object_name, pre_signed_url=pre_signed_url, put_headers=put_headers
        )

    def create_multipart_upload(
        self, bucket_name: str, object_name: str, content_length: int, part_size: int = S3_MULTIPART_CHUNK_SIZE
    ) -> CreateMultipartUploadResponse:
        """
        Tell S3 to start a multipart upload.
        S3 gives us back an upload id, which is included in each part we then upload.

        This action is performed by service account so we use the s3_client not the presigning client.
        """
        response = self.s3_client.create_multipart_upload(Bucket=bucket_name, Key=object_name)
        number_of_parts = ceil(content_length / part_size)
        return CreateMultipartUploadResponse(
            name=object_name,
            upload_id=response["UploadId"],
            number_of_parts=number_of_parts,
            part_size=part_size,
        )

    def create_presigned_upload_part_urls(
        self,
        bucket_name: str,
        object_name: str,
        upload_id: str,
        parts_range_start: int,
        parts_range_end: int,
        md5_checksums: list[str] | None,
    ) -> list[PresignedUploadPartUrlResponse]:
        """Generates a list of pre-signed URLs for a range of parts."""
        urls = []
        for i, part_number in enumerate(range(parts_range_start, parts_range_end + 1)):
            md5_hash = md5_checksums[i] if md5_checksums else None
            params = {
                "Bucket": bucket_name,
                "Key": object_name,
                "UploadId": upload_id,
                "PartNumber": part_number,
            }
            headers = {}
            if md5_hash:
                params["ContentMD5"] = md5_hash
                headers["Content-MD5"] = md5_hash  # used by client

            url = self.s3_pre_signing_client.generate_presigned_url(
                ClientMethod="upload_part",
                ExpiresIn=MULTI_PART_UPLOAD_URL_EXPIRATION_SECONDS,
                Params=params,
            )
            urls.append(PresignedUploadPartUrlResponse(part_number=part_number, pre_signed_url=url, headers=headers))
        return urls

    def complete_multipart_upload(
        self, bucket_name: str, object_name: str, upload_id: str, parts: list[UploadedPart]
    ) -> CompleteMultipartUploadResponse:
        """
        Complete a multipart upload by telling S3 to assemble all previously uploaded parts.
        Each part must have been successfully uploaded before calling this.
        """
        s3_formatted_parts = []
        for part in parts:
            s3_formatted_parts.append({"ETag": part.etag, "PartNumber": part.part_number})
        # Parts have to provided in ascending order (client can upload in whatever order).
        s3_formatted_parts.sort(key=lambda x: x["PartNumber"])

        response = self.s3_client.complete_multipart_upload(
            Bucket=bucket_name,
            Key=object_name,
            UploadId=upload_id,
            MultipartUpload={"Parts": s3_formatted_parts},
        )

        # ETag is the combined etag for the entire object,
        # If differs from the individual part etags and has a suffix like "-N" where N is number of parts.
        # https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/s3/client/complete_multipart_upload.html
        return CompleteMultipartUploadResponse(
            name=object_name,
            version_id=response["VersionId"],
            md5_hash=response["ETag"],
        )

    def abort_multipart_upload(
        self, bucket_name: str, object_name: str, upload_id: str
    ) -> AbortMultipartUploadResponse:
        """
        Abort a multipart upload, deleting any parts that have already been uploaded.
        """
        try:
            self.s3_client.abort_multipart_upload(
                Bucket=bucket_name,
                Key=object_name,
                UploadId=upload_id,
            )
        except self.s3_client.exceptions.NoSuchUpload:
            logger.warning(
                f"Attempted to abort non-existent multipart upload to bucket: {bucket_name}, "
                f"object={object_name}, upload_id={upload_id}\n"
                "Continuing without error."
            )
        return AbortMultipartUploadResponse(name=object_name, upload_id=upload_id)


@lru_cache()
def get_pre_signed_service() -> S3PreSignedService:
    """Dependency to get pre-signed S3 service. To be used in FastAPI endpoints."""
    return S3PreSignedService()
