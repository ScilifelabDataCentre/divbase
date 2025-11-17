"""
S3 files routes. Provides a mix of direct responses and pre-signed URLs depending on the operation.

NOTE:
- The project_name needs to be provided to all routes, it is used by the dependency get_project_member
- Each route can assume the user exists and has access to the project, BUT
we need to use has_required_role to check if they have permission to do the operation.

TODO:
Could be nice to have a detailed list route (so version IDs, sizes, last modified etc).
"""

import logging
from typing import Annotated

from fastapi import APIRouter, Depends, status

from divbase_api.api_config import settings
from divbase_api.crud.projects import has_required_role
from divbase_api.crud.s3 import check_files_already_exist_by_checksum
from divbase_api.deps import get_project_member
from divbase_api.exceptions import AuthorizationError, TooManyObjectsInRequestError
from divbase_api.models.projects import ProjectDB, ProjectRoles
from divbase_api.models.users import UserDB
from divbase_api.services.pre_signed_urls import S3PreSignedService, get_pre_signed_service
from divbase_api.services.s3_client import S3FileManager
from divbase_lib.api_schemas.s3 import (
    CheckFileExistsRequest,
    DownloadObjectRequest,
    ExistingFileResponse,
    PreSignedDownloadResponse,
    PreSignedUploadResponse,
    UploadObjectRequest,
)

logger = logging.getLogger(__name__)

s3_router = APIRouter()


def check_too_many_objects_in_request(numb_objects: int, max_objects: int = 100) -> None:
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

    check_too_many_objects_in_request(len(objects_to_download))

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

    return s3_file_manager.list_files(bucket_name=project.bucket_name)


@s3_router.post("/upload", status_code=status.HTTP_200_OK, response_model=list[PreSignedUploadResponse])
async def generate_upload_url(
    project_name: str,
    objects: list[UploadObjectRequest],
    s3_signer_service: Annotated[S3PreSignedService, Depends(get_pre_signed_service)],
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
):
    """Generate pre-signed POST urls to upload 1 or more file to S3. Max 100 files at a time."""
    project, current_user, role = project_and_user_and_role
    if not has_required_role(role, ProjectRoles.EDIT):
        raise AuthorizationError("You don't have permission to upload files to this project.")

    check_too_many_objects_in_request(len(objects))

    response = []
    for obj in objects:
        pre_signed_response = s3_signer_service.create_presigned_url_for_upload(
            bucket_name=project.bucket_name,
            object_name=obj.name,
            md5_hash=obj.md5_hash,
        )
        response.append(pre_signed_response)

    return response


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

    check_too_many_objects_in_request(len(objects))

    s3_file_manager = S3FileManager(
        url=settings.s3.endpoint_url,
        access_key=settings.s3.access_key.get_secret_value(),
        secret_key=settings.s3.secret_key.get_secret_value(),
    )

    return s3_file_manager.soft_delete_objects(objects=objects, bucket_name=project.bucket_name)


# using POST as GET with body is not considered good practice
@s3_router.post("/check-exists", status_code=status.HTTP_200_OK, response_model=list[ExistingFileResponse])
async def check_file_already_exists_by_checksum(
    files_to_check: list[CheckFileExistsRequest],
    project_name: str,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
):
    """
    Check if the files provided already exist in the project's bucket.
    Compares the MD5 checksums of files of the same name. Returns those that do.
    Max 100 files at a time.
    """
    project, current_user, role = project_and_user_and_role
    if not has_required_role(role, ProjectRoles.READ):
        raise AuthorizationError("You don't have permission to check files in this project.")

    check_too_many_objects_in_request(len(files_to_check))

    return check_files_already_exist_by_checksum(files_to_check=files_to_check, bucket_name=project.bucket_name)
