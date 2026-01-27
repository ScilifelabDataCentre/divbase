"""
Crud operations for the S3 endpoints
"""

from divbase_api.api_config import settings
from divbase_api.services.s3_client import S3FileManager
from divbase_lib.api_schemas.s3 import FileChecksumResponse


def get_s3_checksums(bucket_name: str, files_to_check: list[str]) -> list[FileChecksumResponse]:
    """
    Given a list of object names, return a list of their names and checksums for those that exist in S3.

    TODO: If we have a bucket with a very large number of files, this will result in a lot of calls to S3.
    But this is batched at blocks of 100 in the route layer, so should be manageable for now.
    If we consider a list_objects_v2 approach we have to think carefully about pagination there,
    other we make a lot of calls to S3 anyway.
    """
    s3_file_manager = S3FileManager(
        url=settings.s3.endpoint_url,
        access_key=settings.s3.access_key.get_secret_value(),
        secret_key=settings.s3.secret_key.get_secret_value(),
    )

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
