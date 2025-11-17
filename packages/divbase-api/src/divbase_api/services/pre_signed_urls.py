"""
Logic to generate presigned URLs for S3 objects.
Gives end users access to S3 objects without needing AWS credentials.

Certain operations like "list files" etc... instead performed directly by the backend and returned to end user.
Pre-signed url approach not used when the request to s3 can be done in the (user to API) request-response cycle.

TODO, what is user wants to download same object at multiple versions? - will approach fail as dict keys must be unique.
"""

import logging
from functools import lru_cache

import boto3

from divbase_api.api_config import settings
from divbase_lib.api_schemas.s3 import PreSignedDownloadResponse, PreSignedUploadResponse

logger = logging.getLogger(__name__)


class S3PreSignedService:
    """
    Service to create pre-signed URLs for S3 object upload/download.
    Knows nothing about user authentication/authorization.
    """

    def __init__(self):
        self.s3_client = boto3.client(
            "s3",
            endpoint_url=settings.s3.presigning_url,
            aws_access_key_id=settings.s3.access_key.get_secret_value(),
            aws_secret_access_key=settings.s3.secret_key.get_secret_value(),
        )

    def create_presigned_url_for_download(
        self, bucket_name: str, object_name: str, version_id: str | None
    ) -> PreSignedDownloadResponse:
        """
        Generate a presigned URL for S3 object download

        The generate_presigned_url method from boto3 has following params:
            :param client_method_name: Name of the S3.Client method, e.g., 'list_buckets'
            :param method_parameters: Dictionary of parameters to send to the method
            :param expiration: Time in seconds for the presigned URL to remain valid
            :return: Presigned URL as string. If error, returns None.
        """
        extra_args = {"VersionId": version_id} if version_id else None
        url = self.s3_client.generate_presigned_url(
            ClientMethod="get_object",
            Params={"Bucket": bucket_name, "Key": object_name, **(extra_args or {})},
            ExpiresIn=3600,  # 1 hour
        )

        return PreSignedDownloadResponse(
            name=object_name,
            pre_signed_url=url,
            version_id=version_id,
        )

    def create_presigned_url_for_upload(
        self, bucket_name: str, object_name: str, content_length: int, md5_hash: str | None = None
    ) -> PreSignedUploadResponse:
        """
        Generate a presigned S3 PUT URL to upload a file to S3.
        The reponse object contains the object name, pre-signed URL, and any headers that must be included in the PUT request.

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

        pre_signed_url = self.s3_client.generate_presigned_url(
            HttpMethod="PUT",
            ClientMethod="put_object",
            Params=upload_args,
            ExpiresIn=3600 * 24,  # 24 hours
        )
        return PreSignedUploadResponse(name=object_name, pre_signed_url=pre_signed_url, put_headers=put_headers)


@lru_cache()
def get_pre_signed_service() -> S3PreSignedService:
    """Dependency to get pre-signed S3 service. To be used in FastAPI endpoints."""
    return S3PreSignedService()
