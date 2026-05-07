"""
Crud operations for the S3 endpoints
"""

import logging

import boto3
from botocore.exceptions import ClientError
from pydantic import SecretStr

from divbase_api.api_config import api_settings
from divbase_api.services.pre_signed_urls import S3PreSignedService
from divbase_api.services.s3_client import S3FileManager
from divbase_lib.api_schemas.s3 import FileChecksumResponse

logger = logging.getLogger(__name__)


def get_s3_file_manager() -> S3FileManager:
    """Dependency to get an S3FileManager instance. Dependency injected into FastAPI endpoints."""
    return S3FileManager(
        url=api_settings.s3.endpoint_url,
        access_key=api_settings.s3.access_key.get_secret_value(),
        secret_key=api_settings.s3.secret_key.get_secret_value(),
    )


def get_pre_signed_service() -> S3PreSignedService:
    """Dependency to get pre-signed S3 service. Dependency injected into FastAPI endpoints."""
    return S3PreSignedService()


def get_s3_checksums(bucket_name: str, files_to_check: list[str]) -> list[FileChecksumResponse]:
    """
    Given a list of object names, return a list of their names and checksums for those that exist in S3.

    NOTE: If we have a bucket with a very large number of files, this will result in a lot of calls to S3.
    But this is batched at blocks of 100 in the route layer, so should be manageable for now.
    If we consider a list_objects_v2 approach we have to think carefully about pagination there,
    other we make a lot of calls to S3 anyway.
    """
    s3_file_manager = get_s3_file_manager()

    response = []
    for file_name in files_to_check:
        s3_checksum = s3_file_manager.get_object_checksum_if_exists(bucket_name=bucket_name, object_name=file_name)
        if s3_checksum:
            response.append(
                FileChecksumResponse(
                    object_name=file_name,
                    md5_checksum=s3_checksum,
                )
            )

    return response


def validate_s3_service_account(
    endpoint_url: str, bucket_prefix: str, access_key: SecretStr, secret_key: SecretStr
) -> None:
    """
    Validate the S3 service account can connect and has at least some of its expected permissions, raises an error if not.
    Called at API startup (lifespan event) and celery worker startup

    The bucket_prefix defines the prefix of the buckets the service account is supposed to be able to access in each environment:
    for e.g. local dev and test it is:
    divbase-local-{*} , where * is a number starting from 1.
    """
    s3_client = boto3.client(
        "s3",
        endpoint_url=endpoint_url,
        aws_access_key_id=access_key.get_secret_value(),
        aws_secret_access_key=secret_key.get_secret_value(),
    )
    # does not need to exist as permissions checks will happen before a 404
    made_up_bucket_name = f"{bucket_prefix}1s2112asa1231"

    try:
        _ = s3_client.head_object(
            Bucket=made_up_bucket_name,
            Key="a-non-existent-key-for-connection-and-permissions-check",
        )
    except ClientError as e:
        error_code = e.response["Error"]["Code"]
        if error_code == "403":
            raise RuntimeError(
                "S3 service account does not have the correct permissions or access key and secret key."
            ) from e
        if error_code != "404":
            logger.error(
                f"s3.head_object failed but not due to permissions or not found error. Error was: {error_code}"
            )
            raise RuntimeError(
                "S3 service account cannot connect to S3 or seems to not have expected permissions."
            ) from e

    try:
        s3_client.delete_object(
            Bucket=made_up_bucket_name,
            Key="a-non-existent-key-for-connection-and-permissions-check",
            VersionId="00000000000000000000000000000000",
        )
    except ClientError as e:
        error_code = e.response["Error"]["Code"]
        if error_code == "AccessDenied":
            # we expect the service account to not have hard delete permissions s3:DeleteObjectVersion
            return None
        else:
            raise RuntimeError(
                "Permissions check failed: S3 service MAY have s3:DeleteObjectVersion (hard delete permissions)."
            ) from e
    raise RuntimeError(
        "Permissions check failed: S3 service account seems to have s3:DeleteObjectVersion (hard delete permissions)."
    )
