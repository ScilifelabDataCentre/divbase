"""
Routes for users to manage their project versioning.

Project versions are the state of all files in a project at a given time point.
These are user defined points in time (like a checkpoint/commit) that the user can refer back to later.
"""

import logging

from fastapi import APIRouter, Depends, status
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.api_config import settings
from divbase_api.crud.project_versions import (
    add_project_version,
    get_project_version_details,
    list_project_versions,
    soft_delete_version,
)
from divbase_api.crud.projects import has_required_role
from divbase_api.db import get_db
from divbase_api.deps import get_project_member
from divbase_api.exceptions import AuthorizationError
from divbase_api.models.projects import ProjectDB, ProjectRoles
from divbase_api.models.users import UserDB
from divbase_api.services.s3_client import S3FileManager
from divbase_lib.api_schemas.project_versions import (
    AddVersionRequest,
    AddVersionResponse,
    DeleteVersionRequest,
    DeleteVersionResponse,
    ProjectVersionDetailResponse,
    ProjectVersionInfo,
)

logger = logging.getLogger(__name__)

project_version_router = APIRouter()


@project_version_router.patch("/add", status_code=status.HTTP_200_OK, response_model=AddVersionResponse)
async def add_version_endpoint(
    project_name: str,
    version_request: AddVersionRequest,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
    db: AsyncSession = Depends(get_db),
):
    """
    Add a new entry to the project versioning db table.

    The entry specifies the current state of all files in the project.
    """
    project, current_user, role = project_and_user_and_role
    if not has_required_role(role, ProjectRoles.EDIT):
        raise AuthorizationError("You don't have permission to add a new project version to this project.")

    s3_file_manager = S3FileManager(
        url=settings.s3.endpoint_url,
        access_key=settings.s3.access_key.get_secret_value(),
        secret_key=settings.s3.secret_key.get_secret_value(),
    )

    new_version = await add_project_version(
        db=db,
        name=version_request.name,
        description=version_request.description,
        project=project,
        user_id=current_user.id,
        s3_file_manager=s3_file_manager,
    )

    return new_version


@project_version_router.get("/list", status_code=status.HTTP_200_OK, response_model=list[ProjectVersionInfo])
async def list_versions_endpoint(
    project_name: str,
    include_deleted: bool = False,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
    db: AsyncSession = Depends(get_db),
):
    """
    List all project versions for the project.
    Returns basic info about each version (name, description, timestamp, is_deleted),
    but not a list of all files at that version.

    You can also return soft-deleted versions by including the flag 'include_deleted'.
    """
    project, current_user, role = project_and_user_and_role
    if not has_required_role(role, ProjectRoles.READ):
        raise AuthorizationError("You don't have permission to list project versions for this project.")

    return await list_project_versions(db=db, project_id=project.id, include_deleted=include_deleted)


@project_version_router.get(
    "/version_details", status_code=status.HTTP_200_OK, response_model=ProjectVersionDetailResponse
)
async def project_version_details_endpoint(
    project_name: str,
    version_name: str,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
    db: AsyncSession = Depends(get_db),
):
    """
    List all files and their version IDs at a specific version of the project.
    """
    project, current_user, role = project_and_user_and_role
    if not has_required_role(role, ProjectRoles.READ):
        raise AuthorizationError("You don't have permission to list files at a project version for this project.")

    return await get_project_version_details(
        db=db,
        project_id=project.id,
        version_name=version_name,
    )


@project_version_router.delete("/delete", status_code=status.HTTP_200_OK, response_model=DeleteVersionResponse)
async def delete_version_endpoint(
    project_name: str,
    delete_version_request: DeleteVersionRequest,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
    db: AsyncSession = Depends(get_db),
):
    """
    Delete a version entry in the project versioning table.
    This does not delete the files themselves.
    Deleted version entries older than 30 days will be permanently deleted.
    You can ask a DivBase admin to restore a deleted version within that time period.
    """
    project, current_user, role = project_and_user_and_role
    if not has_required_role(role, ProjectRoles.EDIT):
        raise AuthorizationError("You don't have permission to delete a project version for this project.")
    return await soft_delete_version(
        db=db,
        project_id=project.id,
        version_name=delete_version_request.version_name,
    )
