"""
API routes for task history operations.
"""

import logging

from fastapi import APIRouter, Depends, status
from sqlalchemy.ext.asyncio import AsyncSession
from typing_extensions import Annotated

from divbase_api.crud.projects import has_required_role
from divbase_api.db import get_db
from divbase_api.deps import get_current_user, get_project_member
from divbase_api.exceptions import AuthorizationError
from divbase_api.models.projects import ProjectDB, ProjectRoles
from divbase_api.models.users import UserDB
from divbase_api.services.task_history import (
    get_project_task_history,
    get_task_history_by_id,
    get_user_and_project_task_history,
    get_user_task_history,
)
from divbase_lib.schemas.task_history import TaskHistoryResults

logger = logging.getLogger(__name__)

task_history_router = APIRouter()


@task_history_router.get("/tasks/user", status_code=status.HTTP_200_OK, response_model=TaskHistoryResults)
async def get_all_tasks_for_user(
    current_user: Annotated[UserDB, Depends(get_current_user)],
    db: AsyncSession = Depends(get_db),
) -> TaskHistoryResults:
    """
    Get the task history for the current user. Admin users can view all tasks, non-admin users can only view their own tasks.
    """
    result = await get_user_task_history(
        db=db,
        user_id=current_user.id,
        is_admin=current_user.is_admin,
    )

    result.user_email = current_user.email
    return result


@task_history_router.get(
    "/tasks/user/projects/{project_name}", status_code=status.HTTP_200_OK, response_model=TaskHistoryResults
)
async def get_all_tasks_for_user_and_project(
    project_name: str,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
    db: AsyncSession = Depends(get_db),
) -> TaskHistoryResults:
    """
    Get the task history for the current user and project. Admin users can view all tasks of the project, non-admin users can only view their own tasks of the project.
    """

    project, current_user, role = project_and_user_and_role

    # TODO harmonize how to handle when user has not submitted any tasks, or do not have permissions. A Read user cannot submit tasks...
    if not has_required_role(role, ProjectRoles.EDIT):
        raise AuthorizationError(
            "Project not found or you don't have permission to view task history from this project."
        )

    result = await get_user_and_project_task_history(
        db=db,
        project_id=project.id,
        user_id=current_user.id,
        is_admin=current_user.is_admin,
    )

    result.user_email = current_user.email
    return result


@task_history_router.get("/tasks/{task_id}", status_code=status.HTTP_200_OK, response_model=TaskHistoryResults)
async def get_task_by_id(
    task_id: str,
    current_user: Annotated[UserDB, Depends(get_current_user)],
    db: AsyncSession = Depends(get_db),
) -> TaskHistoryResults:
    """
    Get the history of a specific task ID. Admin users can view any task, non-admin users can only view their own tasks.
    """
    return await get_task_history_by_id(
        db=db,
        task_id=task_id,
        user_id=current_user.id,
        is_admin=current_user.is_admin,
    )


@task_history_router.get("/projects/{project_name}", status_code=status.HTTP_200_OK, response_model=TaskHistoryResults)
async def get_project_tasks(
    project_name: str,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
    db: AsyncSession = Depends(get_db),
) -> TaskHistoryResults:
    """
    Get the task history for a project. Requires MANAGE role or higher.
    """

    project, _, role = project_and_user_and_role

    if not has_required_role(role, ProjectRoles.MANAGE):
        raise AuthorizationError(
            "Project not found or you don't have permission to view task history for this whole project."
        )

    return await get_project_task_history(
        db=db,
        project_id=project.id,
    )
