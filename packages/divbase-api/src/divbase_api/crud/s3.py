"""
Crud operations for the S3 endpoints
"""

from divbase_api.api_config import settings
from divbase_lib.api_schemas.s3 import CheckFileExistsRequest, ExistingFileResponse
from divbase_lib.s3_client import S3FileManager


def check_files_already_exist_by_checksum(
    files_to_check: list[CheckFileExistsRequest], bucket_name: str
) -> list[ExistingFileResponse]:
    """
    Check if files already exist in S3 with the expected MD5 checksums.
    Returns those that do.
    """
    s3_file_manager = S3FileManager(
        url=settings.s3.s3_internal_url,
        access_key=settings.s3.access_key.get_secret_value(),
        secret_key=settings.s3.secret_key.get_secret_value(),
    )

    response = []
    for file in files_to_check:
        # TODO - this wont scale well for a large number of files.p
        # Either run in parallel or use list_objects command but cleverly so we don't pull down everything if lots of files in the bucket.
        checksum = s3_file_manager.get_file_checksum(bucket_name=bucket_name, object_name=file.object_name)
        if checksum and checksum == file.md5_checksum:
            response.append(
                ExistingFileResponse(
                    object_name=file.object_name,
                    md5_checksum=file.md5_checksum,
                    matching_object_name=file.object_name,
                )
            )

    return response
