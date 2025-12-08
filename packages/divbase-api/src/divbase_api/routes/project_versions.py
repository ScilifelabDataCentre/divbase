"""
Routes for users to manage the bucket versioning file in their projects bucket.
TODO - update doc string.

The BucketVersionManager service handles internal logic and send any changes up to S3.
"""

import logging

from fastapi import APIRouter, Depends, status
from fastapi.concurrency import run_in_threadpool

from divbase_api.crud.projects import has_required_role
from divbase_api.deps import get_project_member
from divbase_api.exceptions import AuthorizationError
from divbase_api.models.projects import ProjectDB, ProjectRoles
from divbase_api.models.users import UserDB
from divbase_api.services.bucket_versions import create_bucket_version_manager
from divbase_lib.api_schemas.project_versions import (
    AddVersionRequest,
    AddVersionResponse,
    CreateVersioningFileRequest,
    CreateVersioningFileResponse,
    DeleteVersionRequest,
    DeleteVersionResponse,
    FilesAtVersionResponse,
    VersionListResponse,
)

logger = logging.getLogger(__name__)

project_version_router = APIRouter()


@project_version_router.post(
    "/create", status_code=status.HTTP_201_CREATED, response_model=CreateVersioningFileResponse
)
async def create_versioning_file(
    project_name: str,
    version_creation_request: CreateVersioningFileRequest,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
):
    """
    Create a new bucket versioning file for the project.
    """
    project, current_user, role = project_and_user_and_role
    if not has_required_role(role, ProjectRoles.EDIT):
        raise AuthorizationError("You don't have permission to create a bucket versioning file for this project.")

    bucket_manager = create_bucket_version_manager(bucket_name=project.bucket_name)
    return await run_in_threadpool(
        bucket_manager.create_metadata_file,
        name=version_creation_request.name,
        description=version_creation_request.description,
    )


@project_version_router.patch("/add", status_code=status.HTTP_200_OK, response_model=AddVersionResponse)
async def add_version(
    project_name: str,
    version_request: AddVersionRequest,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
):
    """
    Add an entry to the bucket versioning file.

    Specifying the current state of all files in the project's storage bucket.
    """
    project, current_user, role = project_and_user_and_role
    if not has_required_role(role, ProjectRoles.EDIT):
        raise AuthorizationError("You don't have permission to add a new bucket version to this project.")

    bucket_manager = create_bucket_version_manager(bucket_name=project.bucket_name)
    return bucket_manager.add_version(name=version_request.name, description=version_request.description)


@project_version_router.get("/list", status_code=status.HTTP_200_OK, response_model=VersionListResponse)
async def list_versions(
    project_name: str,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
):
    """
    List all versions in the bucket versioning file.
    """
    project, current_user, role = project_and_user_and_role
    if not has_required_role(role, ProjectRoles.READ):
        raise AuthorizationError("You don't have permission to list bucket versions for this project.")

    bucket_manager = create_bucket_version_manager(bucket_name=project.bucket_name)
    return await run_in_threadpool(bucket_manager.get_version_info)


@project_version_router.get("/list_detailed", status_code=status.HTTP_200_OK, response_model=FilesAtVersionResponse)
async def list_files_at_version(
    project_name: str,
    bucket_version: str,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
):
    """
    List all files and their version IDs at a specific version of the bucket, as defined by the bucket versioning file.
    """
    project, current_user, role = project_and_user_and_role
    if not has_required_role(role, ProjectRoles.READ):
        raise AuthorizationError("You don't have permission to list files at a bucket version for this project.")

    bucket_manager = create_bucket_version_manager(bucket_name=project.bucket_name)
    return await run_in_threadpool(bucket_manager.all_files_at_bucket_version, bucket_version=bucket_version)


@project_version_router.delete("/delete", status_code=status.HTTP_200_OK, response_model=DeleteVersionResponse)
async def delete_version(
    project_name: str,
    delete_version_request: DeleteVersionRequest,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
):
    """
    Delete an entry in the bucket versioning file specfying a specific state of all files in the project's storage bucket.
    Does not delete the files themselves.
    """
    project, current_user, role = project_and_user_and_role
    if not has_required_role(role, ProjectRoles.EDIT):
        raise AuthorizationError("You don't have permission to delete a bucket version for this project.")

    bucket_manager = create_bucket_version_manager(bucket_name=project.bucket_name)
    deleted_version = await run_in_threadpool(
        bucket_manager.delete_version, bucket_version=delete_version_request.version_name
    )

    return DeleteVersionResponse(deleted_version=deleted_version)
