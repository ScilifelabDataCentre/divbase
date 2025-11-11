import logging
from typing import Any

import httpx
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.api_config import settings
from divbase_api.crud.task_history import (
    get_allowed_task_ids_for_project,
    get_allowed_task_ids_for_user,
    get_allowed_task_ids_for_user_and_project,
)
from divbase_lib.schemas.task_history import FlowerTaskResult, TaskHistoryResults

logger = logging.getLogger(__name__)


# TODO ideally, the flower API could be queried with a list of tasks, but that seem not to be the case. for now, set a limit of 500 tasks to not put a lot of overhead on the call
API_LIMIT = min(100, 500)
REQUEST_URL_WITH_LIMIT = f"{settings.flower.url}/api/tasks?limit={API_LIMIT}"

# TODO make a workaround to check if all allowed ids were returned? could call the Flower API task_ID by task_ID... but it would be inefficient


async def get_user_task_history(
    db: AsyncSession,
    user_id: int,
    project_name: str | None = None,
    is_admin: bool = False,
) -> TaskHistoryResults:
    """
    Get a list of the task history from the Flower API.

    For the case of a results backend purge (task not in flower API results):
    allowed_task_ids uses a db lookup, but all_tasks is fetched from the Flower API.
    Thus, if a task is purged in the results backend, it is naturally excluded.
    """

    allowed_task_ids = await get_allowed_task_ids_for_user(db, user_id, is_admin)

    if not allowed_task_ids:
        return TaskHistoryResults(tasks={})

    all_tasks = _make_flower_request(REQUEST_URL_WITH_LIMIT)

    filtered_results = _filter_flower_results_by_allowed_task_ids(
        all_tasks=all_tasks, allowed_task_ids=allowed_task_ids
    )
    return filtered_results


async def get_user_and_project_task_history(
    db: AsyncSession,
    user_id: int,
    project_id: int | None = None,
    is_admin: bool = False,
    display_limit: int = 50,
) -> TaskHistoryResults:
    """
    Get a list of the task history from the Flower API for a user and project.

    For the case of a results backend purge (task not in flower API results):
    allowed_task_ids uses a db lookup, but all_tasks is fetched from the Flower API.
    Thus, if a task is purged in the results backend, it is naturally excluded.
    """
    allowed_task_ids = await get_allowed_task_ids_for_user_and_project(db, user_id, project_id, is_admin)

    if not allowed_task_ids:
        return TaskHistoryResults(tasks={})

    all_tasks = _make_flower_request(REQUEST_URL_WITH_LIMIT)

    filtered_results = _filter_flower_results_by_allowed_task_ids(
        all_tasks=all_tasks, allowed_task_ids=allowed_task_ids
    )
    return filtered_results


async def get_project_task_history(
    db: AsyncSession,
    project_id: int,
) -> TaskHistoryResults:
    """
    Get the the task history of a project from the Flower API.

    """

    allowed_task_ids = await get_allowed_task_ids_for_project(db, project_id)

    if not allowed_task_ids:
        return TaskHistoryResults(tasks={})

    all_tasks = _make_flower_request(REQUEST_URL_WITH_LIMIT)

    filtered_results = _filter_flower_results_by_allowed_task_ids(
        all_tasks=all_tasks, allowed_task_ids=allowed_task_ids
    )
    return filtered_results


async def get_task_history_by_id(
    db: AsyncSession,
    task_id: str,
    user_id: int,
    is_admin: bool = False,
) -> TaskHistoryResults:
    """
    Get the task history from the Flower API for a specific task ID.
    """
    allowed_task_ids = await get_allowed_task_ids_for_user(db, user_id, is_admin)

    if task_id not in allowed_task_ids:
        return TaskHistoryResults(tasks={})

    request_url = f"{settings.flower.url}/api/task/info/{task_id}"
    task_data = _make_flower_request(request_url)

    if not task_data:
        return TaskHistoryResults(tasks={})
    else:
        return TaskHistoryResults(tasks={task_id: FlowerTaskResult(**task_data)})


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


def _filter_flower_results_by_allowed_task_ids(
    all_tasks: dict[str, Any],
    allowed_task_ids: set[str],
) -> TaskHistoryResults:
    """Filter the Flower API results to include only allowed task IDs."""

    filtered_tasks = {}
    for tid, task_data in all_tasks.items():
        if tid in allowed_task_ids:
            filtered_tasks[tid] = FlowerTaskResult(**task_data)

    return TaskHistoryResults(tasks=filtered_tasks)
