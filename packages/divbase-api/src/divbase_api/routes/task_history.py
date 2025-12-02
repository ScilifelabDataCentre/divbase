"""
API routes for task history operations.
"""

import logging

from fastapi import APIRouter, Depends, status
from sqlalchemy.ext.asyncio import AsyncSession
from typing_extensions import Annotated

from divbase_api.crud.projects import check_if_user_is_not_only_read_user_in_all_their_projects, has_required_role
from divbase_api.crud.task_history import get_tasks_pg
from divbase_api.db import get_db
from divbase_api.deps import get_current_user, get_project_member
from divbase_api.exceptions import AuthorizationError
from divbase_api.models.projects import ProjectDB, ProjectRoles
from divbase_api.models.users import UserDB
from divbase_api.services.task_history import (
    deserialize_tasks_to_result,
)
from divbase_lib.api_schemas.task_history import (
    TaskHistoryResults,
)

logger = logging.getLogger(__name__)

task_history_router = APIRouter()

READ_USER_ERROR_MSG = "You do not have access view to task history from any projects. You need to have at least one project where you have an EDIT role or higher."


@task_history_router.get("/tasks/user", status_code=status.HTTP_200_OK, response_model=TaskHistoryResults)
async def get_all_tasks_for_user(
    current_user: Annotated[UserDB, Depends(get_current_user)],
    db: AsyncSession = Depends(get_db),
) -> TaskHistoryResults:
    """
    Get the task history for the current user. Admin users can view all tasks (even if not member of the projects), non-admin users can only view their own tasks.
    """

    serialized_tasks = await get_tasks_pg(
        db=db,
        user_id=current_user.id,
        is_admin=current_user.is_admin,
    )
    result = deserialize_tasks_to_result(serialized_tasks)

    if result.tasks == {}:
        user_has_at_least_one_edit_role = await check_if_user_is_not_only_read_user_in_all_their_projects(
            db=db,
            user_id=current_user.id,
        )

        if not user_has_at_least_one_edit_role:
            raise AuthorizationError(READ_USER_ERROR_MSG)

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
    Get the task history for the current user and project. Admin users can view all tasks of the project (even if not member of the projects), non-admin users can only view their own tasks of the project.
    """

    project, current_user, role = project_and_user_and_role

    if not has_required_role(role, ProjectRoles.EDIT):
        raise AuthorizationError(
            "Project not found or you don't have permission to view task history from this project."
        )

    serialized_tasks = await get_tasks_pg(
        db=db,
        user_id=current_user.id,
        project_id=project.id,
        is_admin=current_user.is_admin,
    )
    result = deserialize_tasks_to_result(serialized_tasks)
    result.user_email = current_user.email
    return result


@task_history_router.get("/projects/{project_name}", status_code=status.HTTP_200_OK, response_model=TaskHistoryResults)
async def get_project_tasks(
    project_name: str,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member),
    db: AsyncSession = Depends(get_db),
) -> TaskHistoryResults:
    """
    Get the task history for a project. Requires MANAGE role or higher. Admin users can view all tasks of the project (even if not member of the projects).
    """

    project, current_user, role = project_and_user_and_role

    if not has_required_role(role, ProjectRoles.MANAGE) and not current_user.is_admin:
        raise AuthorizationError(
            "Project not found or you don't have permission to view task history for this whole project."
        )

    serialized_tasks = await get_tasks_pg(db=db, project_id=project.id)
    return deserialize_tasks_to_result(serialized_tasks)


@task_history_router.get("/tasks/{task_id}", status_code=status.HTTP_200_OK, response_model=TaskHistoryResults)
async def get_task_by_id(
    task_id: str,
    current_user: Annotated[UserDB, Depends(get_current_user)],
    db: AsyncSession = Depends(get_db),
) -> TaskHistoryResults:
    """
    Get the history of a specific task ID.
    - Admin users can view any task.
    - Non-admin users can view tasks they submitted OR tasks from projects where they have MANAGE role.

    Since the request body does not include project ID, permissions checks are made in get_tasks_pg.
    """

    serialized_task = await get_tasks_pg(
        db=db,
        task_id=task_id,
        user_id=current_user.id,
        is_admin=current_user.is_admin,
        require_manager_role=True,
    )

    if not serialized_task:
        raise AuthorizationError("Task ID not found or you don't have permission to view the history for this task ID.")

    result = deserialize_tasks_to_result([serialized_task])
    return result
