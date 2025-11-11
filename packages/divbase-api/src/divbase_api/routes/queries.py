"""
The API server for DivBase.
"""

import logging
import sys

from fastapi import APIRouter, Depends, status
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.api_config import settings
from divbase_api.crud.projects import has_required_role
from divbase_api.crud.task_history import record_pending_task
from divbase_api.db import get_db
from divbase_api.deps import get_project_member
from divbase_api.exceptions import AuthorizationError
from divbase_api.models.projects import ProjectDB, ProjectRoles
from divbase_api.models.users import UserDB
from divbase_api.worker.tasks import (
    bcftools_pipe_task,
    sample_metadata_query_task,
)
from divbase_lib.schemas.queries import (
    BcftoolsQueryKwargs,
    BcftoolsQueryRequest,
    SampleMetadataQueryKwargs,
    SampleMetadataQueryRequest,
)

logging.basicConfig(level=settings.api.log_level, handlers=[logging.StreamHandler(sys.stdout)])

logger = logging.getLogger(__name__)

query_router = APIRouter()

# TODO harmonize function names
# TODO response_model like in the versioning routes


@query_router.post("/sample-metadata/projects/{project_name}", status_code=status.HTTP_200_OK)
async def sample_metadata_query(
    sample_metadata_query_request: SampleMetadataQueryRequest,
    project_name: str,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
    db: AsyncSession = Depends(get_db),
):
    """
    Submit a sample metadata query for the specified project.
    """
    project, current_user, role = project_and_user_and_role

    if not has_required_role(role, ProjectRoles.EDIT):
        raise AuthorizationError("You don't have permission to query this project.")

    task_kwargs = SampleMetadataQueryKwargs(
        tsv_filter=sample_metadata_query_request.tsv_filter,
        metadata_tsv_name=sample_metadata_query_request.metadata_tsv_name,
        bucket_name=project.bucket_name,
        project_id=project.id,
        user_name=current_user.email,
    )

    results = sample_metadata_query_task.apply_async(kwargs=task_kwargs.model_dump())
    await record_pending_task(db=db, task_id=results.id, user_id=current_user.id, project_id=project.id)

    result_dict = results.get(timeout=10)  # TODO think about what happens if this timeout is reached

    if "error" in result_dict:
        error_type = result_dict.get("type", "ServerError")
        error_details = result_dict.get("error", "Unknown error occurred")
        return {"detail": error_details, "type": error_type}

    return {
        "sample_and_filename_subset": result_dict["sample_and_filename_subset"],
        "unique_sample_ids": result_dict["unique_sample_ids"],
        "unique_filenames": result_dict["unique_filenames"],
        "query_message": result_dict["query_message"],
    }


@query_router.post("/bcftools-pipe/projects/{project_name}", status_code=status.HTTP_201_CREATED)
async def create_bcftools_jobs(
    bcftools_query_request: BcftoolsQueryRequest,
    project_name: str,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
    db: AsyncSession = Depends(get_db),
):
    """
    Create a new bcftools query job for the specified project.
    """
    project, current_user, role = project_and_user_and_role

    if not has_required_role(role, ProjectRoles.EDIT):
        raise AuthorizationError("You don't have permission to query this project.")

    task_kwargs = BcftoolsQueryKwargs(
        tsv_filter=bcftools_query_request.tsv_filter,
        command=bcftools_query_request.command,
        metadata_tsv_name=bcftools_query_request.metadata_tsv_name,
        bucket_name=project.bucket_name,
        project_id=project.id,
        user_name=current_user.email,
    )

    results = bcftools_pipe_task.apply_async(kwargs=task_kwargs.model_dump())
    await record_pending_task(db=db, task_id=results.id, user_id=current_user.id, project_id=project.id)
    return results.id
