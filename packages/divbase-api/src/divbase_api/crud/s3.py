"""
Crud operations for the S3 endpoints
"""

from divbase_api.api_config import settings
from divbase_api.services.s3_client import S3FileManager
from divbase_lib.api_schemas.s3 import CheckFileExistsRequest, ExistingFileResponse


def check_files_already_exist_by_checksum(
    files_to_check: list[CheckFileExistsRequest], bucket_name: str
) -> list[ExistingFileResponse]:
    """
    Check if files already exist in S3 with the expected MD5 checksums.
    Returns only those that do.

    TODO:
    The function tries to optimize the number of S3 calls made by using
    a bulk listing approach when there are 5 or more files to check.
    This probably needs some testing as the performance trade-offs will depend
    on the number of files in the bucket vs those being checked and how relatively slow the s3 operations are.
    """
    s3_file_manager = S3FileManager(
        url=settings.s3.s3_endpoint_url,
        access_key=settings.s3.access_key.get_secret_value(),
        secret_key=settings.s3.secret_key.get_secret_value(),
    )

    response: list[ExistingFileResponse] = []
    if len(files_to_check) < 5:
        for file in files_to_check:
            s3_checksum = s3_file_manager.get_file_checksum(bucket_name=bucket_name, object_name=file.object_name)
            if s3_checksum and s3_checksum == file.md5_checksum:
                response.append(
                    ExistingFileResponse(
                        object_name=file.object_name,
                        md5_checksum=file.md5_checksum,
                        matching_object_name=file.object_name,
                    )
                )
    else:
        object_names = [file.object_name for file in files_to_check]
        existing_checksums = s3_file_manager.get_multiple_checksums(bucket_name=bucket_name, object_names=object_names)

        for file in files_to_check:
            s3_checksum = existing_checksums.get(file.object_name)
            if s3_checksum and s3_checksum == file.md5_checksum:
                response.append(
                    ExistingFileResponse(
                        object_name=file.object_name,
                        md5_checksum=file.md5_checksum,
                        matching_object_name=file.object_name,
                    )
                )

    return response
