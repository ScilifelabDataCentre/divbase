"""
API routes for query operations.
"""

import logging
import sys

import celery
from fastapi import APIRouter, Depends, HTTPException, status
from fastapi.concurrency import run_in_threadpool
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.api_config import settings
from divbase_api.crud.projects import has_required_role
from divbase_api.crud.task_history import create_task_history_entry, update_task_history_entry_with_celery_task_id
from divbase_api.db import get_db
from divbase_api.deps import get_project_member
from divbase_api.exceptions import AuthorizationError, VCFDimensionsEntryMissingError
from divbase_api.models.projects import ProjectDB, ProjectRoles
from divbase_api.models.users import UserDB
from divbase_api.worker.tasks import (
    bcftools_pipe_task,
    sample_metadata_query_task,
)
from divbase_lib.api_schemas.queries import (
    BcftoolsQueryKwargs,
    BcftoolsQueryRequest,
    SampleMetadataQueryKwargs,
    SampleMetadataQueryRequest,
    SampleMetadataQueryTaskResult,
)
from divbase_lib.exceptions import (
    SidecarColumnNotFoundError,
    SidecarInvalidFilterError,
    SidecarSampleIDError,
)

logging.basicConfig(level=settings.api.log_level, handlers=[logging.StreamHandler(sys.stderr)])

logger = logging.getLogger(__name__)

query_router = APIRouter()

# TODO harmonize function names


@query_router.post(
    "/sample-metadata/projects/{project_name}",
    status_code=status.HTTP_200_OK,
    response_model=SampleMetadataQueryTaskResult,
)
async def sample_metadata_query(
    sample_metadata_query_request: SampleMetadataQueryRequest,
    project_name: str,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
    db: AsyncSession = Depends(get_db),
) -> SampleMetadataQueryTaskResult:
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
        project_name=project.name,
        user_id=current_user.id,
    )

    results = sample_metadata_query_task.apply_async(kwargs=task_kwargs.model_dump())

    _ = await create_task_history_entry(
        user_id=current_user.id,
        project_id=project.id,
        task_id=results.id,
        db=db,
    )

    try:
        result_dict = await run_in_threadpool(results.get, timeout=10)
        # TODO - consider if we split this into 2 routes to handle time out issues on CLI side.
        # Route 1, create job and get back job id.
        # Route 2, get job result by id (with status etc), CLI can poll until done.

    except (SidecarInvalidFilterError, SidecarColumnNotFoundError, SidecarSampleIDError) as e:
        # Catch validation errors (mixed types, missing columns, invalid Sample_IDs) and return 400
        error_message = str(e)
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=error_message) from None
    except VCFDimensionsEntryMissingError:
        # Catch and raise anew to avoid duplications in the error message
        raise VCFDimensionsEntryMissingError(project_name=project.name) from None
    except celery.exceptions.TimeoutError:  # type: ignore
        error_message = (
            f"The query is still being processed and has Task ID: {results.id}. \n"
            f"Please check back later for the results. \n"
            f"To check the status of the query you can use the following command: \n"
            f"divbase-cli task-history id {results.id}"
        )
        raise HTTPException(status_code=status.HTTP_408_REQUEST_TIMEOUT, detail=error_message) from None
    except FileNotFoundError:
        error_message = (
            f"The sample metadata TSV file named: {sample_metadata_query_request.metadata_tsv_name} was not found in your project {project.name} \n"
            "Please make sure to upload it first ('divbase-cli files upload ...') and try again."
        )
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=error_message) from None

    return SampleMetadataQueryTaskResult(**result_dict)


@query_router.post("/bcftools-pipe/projects/{project_name}", status_code=status.HTTP_201_CREATED)
async def create_bcftools_jobs(
    bcftools_query_request: BcftoolsQueryRequest,
    project_name: str,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
    db: AsyncSession = Depends(get_db),
) -> int:
    """
    Create a new bcftools query job for the specified project.
    """
    project, current_user, role = project_and_user_and_role

    if not has_required_role(role, ProjectRoles.EDIT):
        raise AuthorizationError("You don't have permission to query this project.")

    job_id = await create_task_history_entry(
        user_id=current_user.id,
        project_id=project.id,
        db=db,
    )

    task_kwargs = BcftoolsQueryKwargs(
        tsv_filter=bcftools_query_request.tsv_filter,
        command=bcftools_query_request.command,
        metadata_tsv_name=bcftools_query_request.metadata_tsv_name,
        bucket_name=project.bucket_name,
        project_id=project.id,
        project_name=project.name,
        user_id=current_user.id,
        job_id=job_id,
    )

    result = bcftools_pipe_task.apply_async(kwargs=task_kwargs.model_dump())

    await update_task_history_entry_with_celery_task_id(
        job_id=job_id,
        task_id=result.id,
        db=db,
    )

    return job_id
