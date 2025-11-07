"""
API routes for task history operations.
"""

import ast
import logging

from fastapi import APIRouter, Depends
from typing_extensions import Annotated

from divbase_api.deps import get_current_user
from divbase_api.get_task_history import get_task_history
from divbase_api.models.users import UserDB

logger = logging.getLogger(__name__)

task_history_router = APIRouter()


@task_history_router.get("/list")
def get_all_tasks(current_user: Annotated[UserDB, Depends(get_current_user)], limit: int = 10):
    """
    Get the task history for the current user. Admin users can view all tasks, non-admin users can only view their own tasks.
    """
    all_tasks = get_task_history(display_limit=limit)

    if current_user.is_admin:
        return all_tasks

    filtered_tasks = []
    for task_id, task_data, timestamp in all_tasks:
        kwargs = task_data.get("kwargs", "{}")

        try:
            parsed_kwargs = ast.literal_eval(kwargs)
        except ValueError as e:
            logger.warning(f"Could not parse kwargs: {e}")
            return "Unknown"

        submitter = parsed_kwargs.get("user_name")
        if submitter == current_user.email:
            filtered_tasks.append((task_id, task_data, timestamp))

    return filtered_tasks


@task_history_router.get("/list/{task_id}")
def get_task_by_id(
    task_id: str,
    current_user: Annotated[UserDB, Depends(get_current_user)],
):
    """
    Get the history of a specific task ID. Admin users can view any task, non-admin users can only view their own tasks.
    """
    task_items = get_task_history(task_id=task_id)

    if not current_user.is_admin:
        for _, task_data, _ in task_items:
            kwargs = task_data.get("kwargs", {})

            try:
                parsed_kwargs = ast.literal_eval(kwargs)
            except ValueError as e:
                logger.warning(f"Could not parse kwargs: {e}")
                return "Unknown"

            submitter = parsed_kwargs.get("user_name")
            if submitter != current_user.email:
                return []

    return task_items
