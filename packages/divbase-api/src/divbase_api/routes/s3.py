"""
S3 files routes. Provides a mix of direct responses and pre-signed URLs depending on the operation.

NOTE:
- The project_name needs to be provided to all routes, it is used by the dependency get_project_member
- Each route can assume the user exists and has access to the project, BUT
we need to use has_required_role to check if they have permission to do the operation.
- To avoid blocking the event loop when using the S3 client (boto3 is a sync SDK), we run these operations in a threadpool.

TODO:
Could be nice to have a detailed list route (so version IDs, sizes, last modified etc).
"""

import logging
from typing import Annotated

from fastapi import APIRouter, Depends, HTTPException, status
from fastapi.concurrency import run_in_threadpool

from divbase_api.api_config import settings
from divbase_api.crud.projects import has_required_role
from divbase_api.crud.s3 import get_s3_checksums
from divbase_api.deps import get_project_member
from divbase_api.exceptions import AuthorizationError, TooManyObjectsInRequestError
from divbase_api.models.projects import ProjectDB, ProjectRoles
from divbase_api.models.users import UserDB
from divbase_api.services.pre_signed_urls import S3PreSignedService, get_pre_signed_service
from divbase_api.services.s3_client import S3FileManager
from divbase_lib.api_schemas.s3 import (
    AbortMultipartUploadRequest,
    AbortMultipartUploadResponse,
    CompleteMultipartUploadRequest,
    CompleteMultipartUploadResponse,
    CreateMultipartUploadRequest,
    CreateMultipartUploadResponse,
    DownloadObjectRequest,
    FileChecksumResponse,
    GetPresignedPartUrlsRequest,
    PreSignedDownloadResponse,
    PreSignedSinglePartUploadResponse,
    PresignedUploadPartUrlResponse,
    UploadSinglePartObjectRequest,
)
from divbase_lib.divbase_constants import MAX_S3_API_BATCH_SIZE

logger = logging.getLogger(__name__)

s3_router = APIRouter()


def check_too_many_objects_in_request(
    numb_objects: int,
    max_objects: int = MAX_S3_API_BATCH_SIZE,
) -> None:
    """
    Helper function to check if too many objects are being requested to be worked on
    in one single request.
    """
    if numb_objects > max_objects:
        raise TooManyObjectsInRequestError(
            f"Too many objects/files provided to be worked on. You provided: {numb_objects}. Maximum allowed: {max_objects}."
        )


# Post request instead of GET as GET doesn't support/encourage body content.
@s3_router.post("/download", status_code=status.HTTP_200_OK, response_model=list[PreSignedDownloadResponse])
async def generate_download_urls(
    project_name: str,
    objects_to_download: list[DownloadObjectRequest],
    s3_signer_service: Annotated[S3PreSignedService, Depends(get_pre_signed_service)],
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
):
    """Generate pre-signed URLs for downloading files at specific versions from S3."""
    project, current_user, role = project_and_user_and_role
    if not has_required_role(role, ProjectRoles.READ):
        raise AuthorizationError("You don't have permission to download files from this project.")

    check_too_many_objects_in_request(numb_objects=len(objects_to_download))

    response = []
    for obj in objects_to_download:
        pre_signed_response = s3_signer_service.create_presigned_url_for_download(
            bucket_name=project.bucket_name, object_name=obj.name, version_id=obj.version_id
        )
        response.append(pre_signed_response)

    return response


@s3_router.get("/", status_code=status.HTTP_200_OK, response_model=list[str])
async def list_files(
    project_name: str,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
):
    """List all files in the project's bucket"""
    project, current_user, role = project_and_user_and_role
    if not has_required_role(role, ProjectRoles.READ):
        raise AuthorizationError("You don't have permission to list files in this project.")

    s3_file_manager = S3FileManager(
        url=settings.s3.endpoint_url,
        access_key=settings.s3.access_key.get_secret_value(),
        secret_key=settings.s3.secret_key.get_secret_value(),
    )

    return await run_in_threadpool(s3_file_manager.list_files, bucket_name=project.bucket_name)


@s3_router.post(
    "/upload/single-part", status_code=status.HTTP_200_OK, response_model=list[PreSignedSinglePartUploadResponse]
)
async def generate_single_part_upload_urls(
    project_name: str,
    objects: list[UploadSinglePartObjectRequest],
    s3_signer_service: Annotated[S3PreSignedService, Depends(get_pre_signed_service)],
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
):
    """
    Generate pre-signed POST urls to upload 1 or more file to S3 via single part uploads.
    Constraints:
    - Max 100 files at a time.
    - Each file must be less than 5GB in size.

    Larger files must use multipart upload.
    """
    project, current_user, role = project_and_user_and_role
    if not has_required_role(role, ProjectRoles.EDIT):
        raise AuthorizationError("You don't have permission to upload files to this project.")

    check_too_many_objects_in_request(numb_objects=len(objects))

    response = []
    for obj in objects:
        pre_signed_response = s3_signer_service.create_presigned_url_for_single_part_upload(
            bucket_name=project.bucket_name,
            object_name=obj.name,
            md5_hash=obj.md5_hash,
            content_length=obj.content_length,
        )
        response.append(pre_signed_response)

    return response


@s3_router.post(
    "/upload/multi-part/create", status_code=status.HTTP_200_OK, response_model=CreateMultipartUploadResponse
)
async def create_multi_part_upload(
    project_name: str,
    upload_request: CreateMultipartUploadRequest,
    s3_signer_service: Annotated[S3PreSignedService, Depends(get_pre_signed_service)],
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
):
    project, current_user, role = project_and_user_and_role
    if not has_required_role(role, ProjectRoles.EDIT):
        raise AuthorizationError("You don't have permission to upload files to this project.")

    return await run_in_threadpool(
        s3_signer_service.create_multipart_upload,
        bucket_name=project.bucket_name,
        object_name=upload_request.name,
        content_length=upload_request.content_length,
        part_size=upload_request.part_size,
    )


# using POST as GET with body is not considered good practice
@s3_router.post(
    "/upload/multi-part/part-urls", status_code=status.HTTP_200_OK, response_model=list[PresignedUploadPartUrlResponse]
)
async def get_pre_signed_urls_parts(
    project_name: str,
    parts_request: GetPresignedPartUrlsRequest,
    s3_signer_service: Annotated[S3PreSignedService, Depends(get_pre_signed_service)],
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
):
    project, current_user, role = project_and_user_and_role
    if not has_required_role(role, ProjectRoles.EDIT):
        raise AuthorizationError("You don't have permission to upload files to this project.")

    numb_parts = parts_request.parts_range_end - parts_request.parts_range_start + 1
    check_too_many_objects_in_request(numb_objects=numb_parts)

    if parts_request.md5_checksums and len(parts_request.md5_checksums) != numb_parts:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="The number of md5_checksums must match the number of parts requested.",
        )

    return s3_signer_service.create_presigned_upload_part_urls(
        bucket_name=project.bucket_name,
        object_name=parts_request.name,
        upload_id=parts_request.upload_id,
        parts_range_start=parts_request.parts_range_start,
        parts_range_end=parts_request.parts_range_end,
        md5_checksums=parts_request.md5_checksums,
    )


@s3_router.post(
    "/upload/multi-part/complete", status_code=status.HTTP_200_OK, response_model=CompleteMultipartUploadResponse
)
async def complete_multipart_upload(
    project_name: str,
    complete_request: CompleteMultipartUploadRequest,
    s3_signer_service: Annotated[S3PreSignedService, Depends(get_pre_signed_service)],
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
):
    project, current_user, role = project_and_user_and_role
    if not has_required_role(role, ProjectRoles.EDIT):
        raise AuthorizationError("You don't have permission to upload files to this project.")

    return await run_in_threadpool(
        s3_signer_service.complete_multipart_upload,
        bucket_name=project.bucket_name,
        object_name=complete_request.name,
        upload_id=complete_request.upload_id,
        parts=complete_request.parts,
    )


@s3_router.delete(
    "/upload/multi-part/abort", status_code=status.HTTP_200_OK, response_model=AbortMultipartUploadResponse
)
async def abort_multipart_upload(
    project_name: str,
    abort_request: AbortMultipartUploadRequest,
    s3_signer_service: Annotated[S3PreSignedService, Depends(get_pre_signed_service)],
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
):
    project, current_user, role = project_and_user_and_role
    if not has_required_role(role, ProjectRoles.EDIT):
        raise AuthorizationError("You don't have permission to upload files to this project.")

    s3_signer_service.abort_multipart_upload(
        bucket_name=project.bucket_name,
        object_name=abort_request.name,
        upload_id=abort_request.upload_id,
    )
    return AbortMultipartUploadResponse(name=abort_request.name, upload_id=abort_request.upload_id)


@s3_router.delete("/", status_code=status.HTTP_200_OK, response_model=list[str])
async def soft_delete_files(
    objects: list[str],
    project_name: str,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
):
    """Soft delete files in the project's bucket. This adds a deletion marker to the files, but does not actually delete them."""
    project, current_user, role = project_and_user_and_role
    if not has_required_role(role, ProjectRoles.EDIT):
        raise AuthorizationError("You don't have permission to soft delete files in this project.")

    check_too_many_objects_in_request(numb_objects=len(objects))

    s3_file_manager = S3FileManager(
        url=settings.s3.endpoint_url,
        access_key=settings.s3.access_key.get_secret_value(),
        secret_key=settings.s3.secret_key.get_secret_value(),
    )

    return await run_in_threadpool(
        s3_file_manager.soft_delete_objects, objects=objects, bucket_name=project.bucket_name
    )


# using POST as GET with body is not considered good practice
@s3_router.post("/checksums", status_code=status.HTTP_200_OK, response_model=list[FileChecksumResponse])
async def get_files_checksums(
    object_names: list[str],
    project_name: str,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
):
    """
    Given a list of potential file names in the bucket, return a list of those that exists and their checksums.

    This can be used ahead of time when uploading multiple files and you want to check (before uploading)
    which files already exist in the bucket.

    It is not assumed that all files provided actually exist in the bucket,
    only those that do are returned in the response dict.

    Max 100 files at a time.
    """
    project, current_user, role = project_and_user_and_role
    if not has_required_role(role, ProjectRoles.READ):
        raise AuthorizationError("You don't have permission to get file checksums for this project.")

    check_too_many_objects_in_request(numb_objects=len(object_names))

    return await run_in_threadpool(get_s3_checksums, bucket_name=project.bucket_name, files_to_check=object_names)
