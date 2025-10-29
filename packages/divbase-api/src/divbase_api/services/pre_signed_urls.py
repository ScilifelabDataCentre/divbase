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
from divbase_api.schemas.s3 import PreSignedDownloadResponse, PreSignedUploadResponse

logger = logging.getLogger(__name__)


class S3PreSignedService:
    """
    Service to create pre-signed URLs for S3 object upload/download.
    Knows nothing about user authentication/authorization.
    """

    def __init__(self):
        self.s3_client = boto3.client(
            "s3",
            endpoint_url=settings.s3.s3_external_url,
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
            object_name=object_name,
            pre_signed_url=url,
            version_id=version_id,
        )

    def create_presigned_url_for_upload(
        self, bucket_name: str, object_name: str, fields=None, conditions=None
    ) -> PreSignedUploadResponse:
        """
        Generate a presigned S3 POST URL to upload a file.

        Notes:
        :param fields: Dictionary of prefilled form fields
        :param conditions: List of conditions to include in the policy

        Returns a dictionary with two elements: url and fields. Url is the url to post to.
        Fields is a dictionary filled with the form fields and respective values
        to use when submitting the post.
        """
        response = self.s3_client.generate_presigned_post(
            Bucket=bucket_name,
            Key=object_name,
            Fields=fields,
            Conditions=conditions,
            ExpiresIn=3600 * 24,  # 24 hours
        )
        return PreSignedUploadResponse(
            object_name=object_name,
            post_url=response["url"],
            fields=response["fields"],
        )


@lru_cache()
def get_pre_signed_service() -> S3PreSignedService:
    """Dependency to get pre-signed S3 service. To be used in FastAPI endpoints."""
    return S3PreSignedService()
