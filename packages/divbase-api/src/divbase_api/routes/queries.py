"""
API routes for query operations.
"""

import celery
import structlog
from fastapi import APIRouter, Depends, HTTPException, status
from fastapi.concurrency import run_in_threadpool
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.projects import has_required_role
from divbase_api.crud.queue_status import check_queue_closed_for_new_tasks
from divbase_api.crud.task_history import create_task_history_entry, update_task_history_entry_with_celery_task_id
from divbase_api.crud.vcf_dimensions import get_unique_samples_by_project_async
from divbase_api.db import get_db
from divbase_api.deps import get_project_member
from divbase_api.exceptions import AuthorizationError
from divbase_api.models.projects import ProjectDB, ProjectRoles
from divbase_api.models.users import UserDB
from divbase_api.services.vcf_queries import validate_user_submitted_bcftools_command
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
    DimensionsNotUpToDateWithBucketError,
    SidecarColumnNotFoundError,
    SidecarInvalidFilterError,
    SidecarMetadataFormatError,
    SidecarSampleIDError,
    TaskUserError,
)

logger = structlog.get_logger(__name__)

query_router = APIRouter()

QUERY_AUTHORIZATION_ERROR_MSG = (
    "You don't have permission to query this project. You need at least 'QUERY' level permissions."
)


@query_router.post(
    "/sample-metadata/projects/{project_name}",
    status_code=status.HTTP_200_OK,
    response_model=SampleMetadataQueryTaskResult,
)
async def submit_sample_metadata_query_job_endpoint(
    sample_metadata_query_request: SampleMetadataQueryRequest,
    project_name: str,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
    db: AsyncSession = Depends(get_db),
) -> SampleMetadataQueryTaskResult:
    """
    Submit a sample metadata query for the specified project.
    """
    project, current_user, role = project_and_user_and_role

    if not has_required_role(role, ProjectRoles.QUERY):
        raise AuthorizationError(QUERY_AUTHORIZATION_ERROR_MSG)
    await check_queue_closed_for_new_tasks(db=db, is_admin=current_user.is_admin)

    task_kwargs = SampleMetadataQueryKwargs(
        tsv_filter=sample_metadata_query_request.tsv_filter,
        metadata_tsv_name=sample_metadata_query_request.metadata_tsv_name,
        bucket_name=project.bucket_name,
        project_id=project.id,
        project_name=project.name,
        user_id=current_user.id,
    )

    results = sample_metadata_query_task.apply_async(kwargs=task_kwargs.model_dump())

    job_id = await create_task_history_entry(
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

    except (
        SidecarInvalidFilterError,
        SidecarColumnNotFoundError,
        SidecarSampleIDError,
        SidecarMetadataFormatError,
        TaskUserError,
        DimensionsNotUpToDateWithBucketError,
    ) as e:
        # These are simple exceptions (that inherit from base Exception) that are able to pass through Celery's JSON serialization/deserialization without becoming UnpicklableExceptionWrapper.
        # TaskUserError is a wrapper exception that allow to nest more complex exceptions that would normally trigger UnpicklableExceptionWrapper, and still be able to pass through the serialization/deserialization.
        error_message = str(e)
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=error_message) from None
    except celery.exceptions.TimeoutError:  # type: ignore
        error_message = (
            f"The query is still being processed and has Task ID: {job_id}. \n"
            f"Please check back later for the results. \n"
            f"To check the status of the query you can use the following command: \n"
            f"divbase-cli task-history id {job_id}"
        )
        raise HTTPException(status_code=status.HTTP_408_REQUEST_TIMEOUT, detail=error_message) from None

    return SampleMetadataQueryTaskResult(**result_dict)


@query_router.post("/vcf/projects/{project_name}", status_code=status.HTTP_201_CREATED)
async def submit_vcf_query_job_endpoint(
    bcftools_query_request: BcftoolsQueryRequest,
    project_name: str,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
    db: AsyncSession = Depends(get_db),
) -> int:
    """
    Submit a vcf query job for the specified project.
    """
    project, current_user, role = project_and_user_and_role

    if not has_required_role(role, ProjectRoles.QUERY):
        raise AuthorizationError(QUERY_AUTHORIZATION_ERROR_MSG)
    await check_queue_closed_for_new_tasks(db=db, is_admin=current_user.is_admin)

    if bcftools_query_request.samples is not None:
        project_samples = set(await get_unique_samples_by_project_async(db=db, project_id=project.id))
        mismatched_samples = sorted(set(bcftools_query_request.samples) - project_samples)
        if mismatched_samples:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=(
                    "The following sample IDs were not found in the project's dimensions index: "
                    f"{', '.join(mismatched_samples)}. "
                    "Please verify that the sample names are correctly spelled. "
                    "If the sample names are correct, please make sure the dimensions index is up to date "
                    "by running 'divbase-cli dimensions update --project <project_name>'."
                ),
            )

    try:
        validated_bcftools_commands = validate_user_submitted_bcftools_command(
            command=bcftools_query_request.command,
            all_samples=bcftools_query_request.all_samples,
        )
    except TaskUserError as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e)) from None

    job_id = await create_task_history_entry(
        user_id=current_user.id,
        project_id=project.id,
        db=db,
    )

    task_kwargs = BcftoolsQueryKwargs(
        tsv_filter=bcftools_query_request.tsv_filter,
        command=validated_bcftools_commands,
        metadata_tsv_name=bcftools_query_request.metadata_tsv_name,
        bucket_name=project.bucket_name,
        project_id=project.id,
        project_name=project.name,
        user_id=current_user.id,
        job_id=job_id,
        samples=bcftools_query_request.samples,
        all_samples=bcftools_query_request.all_samples,
    )

    result = bcftools_pipe_task.apply_async(kwargs=task_kwargs.model_dump())

    await update_task_history_entry_with_celery_task_id(
        job_id=job_id,
        task_id=result.id,
        db=db,
    )

    return job_id
