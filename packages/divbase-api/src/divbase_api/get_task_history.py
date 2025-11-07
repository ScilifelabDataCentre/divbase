import logging
from typing import Any

import httpx
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.api_config import settings
from divbase_api.crud.task_history import (
    filter_task_ids_by_project_name,
    get_allowed_task_ids_for_user,
)
from divbase_lib.queries import TaskHistoryResults

logger = logging.getLogger(__name__)


async def get_task_history_list(
    db: AsyncSession,
    user_id: int,
    project_name: str | None = None,
    is_admin: bool = False,
    display_limit: int = 50,
) -> TaskHistoryResults:
    """
    Get a list of the task history from the Flower API.
    """
    api_limit = min(100, display_limit * 5)
    request_url = f"{settings.flower.url}/api/tasks?limit={api_limit}"
    all_tasks = _make_flower_request(request_url)

    allowed_task_ids = await get_allowed_task_ids_for_user(db, user_id, is_admin)

    if project_name:
        allowed_task_ids = await filter_task_ids_by_project_name(db, allowed_task_ids, project_name)

    filtered_tasks = {}
    for tid, data in all_tasks.items():
        if tid in allowed_task_ids:
            filtered_tasks[tid] = data

    return TaskHistoryResults(tasks=filtered_tasks)


async def get_task_history_by_id(
    db: AsyncSession,
    task_id: str,
    user_id: int,
    is_admin: bool = False,
) -> TaskHistoryResults:
    """
    Get the task history from the Flower API for a specific task ID.
    """
    request_url = f"{settings.flower.url}/api/task/info/{task_id}"
    task = _make_flower_request(request_url)

    allowed_task_ids = await get_allowed_task_ids_for_user(db, user_id, is_admin)

    if task_id not in allowed_task_ids:
        return TaskHistoryResults(tasks={})

    return TaskHistoryResults(tasks={task_id: task})


def _make_flower_request(request_url: str) -> dict[str, Any]:
    """
    Make a request to the Flower API for info about tasks.

    Returns the JSON response as a dictionary.
    If multiple tasks returned, you get back a nested dict, where outer keys are task IDs.
    """
    with httpx.Client() as client:
        response = client.get(
            url=request_url,
            timeout=3.0,
            auth=(
                settings.flower.user,
                settings.flower.password.get_secret_value(),
            ),
        )

    if response.status_code != 200:
        raise ConnectionError(f"Failed to fetch tasks info from Flower API. Status code: {response.status_code}")

    return response.json()
