"""
API routes for VCF dimensions (technical metadata) operations.
"""

import logging

from fastapi import APIRouter, Depends, status
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.projects import has_required_role
from divbase_api.crud.task_history import record_pending_task
from divbase_api.crud.vcf_dimensions import (
    get_skipped_vcfs_by_project_async,
    get_vcf_metadata_by_project_async,
)
from divbase_api.db import get_db
from divbase_api.deps import get_project_member
from divbase_api.exceptions import AuthorizationError, VCFDimensionsEntryMissingError
from divbase_api.models.projects import ProjectDB, ProjectRoles
from divbase_api.models.users import UserDB
from divbase_api.worker.tasks import update_vcf_dimensions_task

logger = logging.getLogger(__name__)

vcf_dimensions_router = APIRouter()


@vcf_dimensions_router.get("/projects/{project_name}", status_code=status.HTTP_200_OK)
async def list_vcf_metadata_by_project_name_user_endpoint(
    project_name: str,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
    db: AsyncSession = Depends(get_db),
):
    """
    Get all VCF metadata entries for a project by project name.

    This endpoint is for regular users (CLI/frontend) and uses the standard
    project member authorization pattern.
    """
    project, current_user, role = project_and_user_and_role

    if not has_required_role(role, ProjectRoles.READ):
        raise AuthorizationError("You don't have permission to view VCF dimensions for this project.")

    result = await get_vcf_metadata_by_project_async(db, project.id)
    vcf_files = result.get("vcf_files", [])

    skipped_entries = await get_skipped_vcfs_by_project_async(db, project.id)
    skipped_files = [
        {
            "vcf_file_s3_key": entry.vcf_file_s3_key,
            "s3_version_id": entry.s3_version_id,
            "skip_reason": entry.skip_reason,
        }
        for entry in skipped_entries
    ]

    if not vcf_files and not skipped_files:
        raise VCFDimensionsEntryMissingError(project_name=project.name)

    return {
        "project_id": project.id,
        "project_name": project.name,
        "vcf_file_count": len(vcf_files),
        "vcf_files": vcf_files,
        "skipped_file_count": len(skipped_files),
        "skipped_files": skipped_files,
    }


@vcf_dimensions_router.put("/projects/{project_name}", status_code=status.HTTP_202_ACCEPTED)
async def update_vcf_dimensions_endpoint(
    project_name: str,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
    db: AsyncSession = Depends(get_db),
):
    """
    Update the VCF dimensions files for the specified project
    """
    project, current_user, role = project_and_user_and_role

    if not has_required_role(role, ProjectRoles.EDIT):
        raise AuthorizationError("You don't have permission to update VCF dimensions for this project.")

    task_kwargs = {"bucket_name": project.bucket_name, "project_id": project.id, "user_name": current_user.email}

    results = update_vcf_dimensions_task.apply_async(kwargs=task_kwargs)
    await record_pending_task(db=db, task_id=results.id, user_id=current_user.id, project_id=project.id)
    return results.id
