"""
API routes for task history operations.
"""

import logging

from fastapi import APIRouter, Depends
from typing_extensions import Annotated

from divbase_api.deps import get_current_user
from divbase_api.get_task_history import TaskHistoryResults, get_task_history_by_id, get_task_history_list
from divbase_api.models.users import UserDB

logger = logging.getLogger(__name__)

task_history_router = APIRouter()


@task_history_router.get("/list")
def get_all_tasks(current_user: Annotated[UserDB, Depends(get_current_user)], limit: int) -> TaskHistoryResults:
    """
    Get the task history for the current user. Admin users can view all tasks, non-admin users can only view their own tasks.
    """
    return get_task_history_list(
        display_limit=limit,
        user_id=current_user.id,
        is_admin=current_user.is_admin,
    )


@task_history_router.get("/list/{task_id}")
def get_task_by_id(
    task_id: str,
    current_user: Annotated[UserDB, Depends(get_current_user)],
) -> TaskHistoryResults:
    """
    Get the history of a specific task ID. Admin users can view any task, non-admin users can only view their own tasks.
    """
    return get_task_history_by_id(
        task_id=task_id,
        user_id=current_user.id,
        is_admin=current_user.is_admin,
    )
