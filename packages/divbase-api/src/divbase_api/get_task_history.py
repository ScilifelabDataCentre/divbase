import logging
from typing import Any

import httpx

from divbase_api.api_config import settings
from divbase_api.models.task_history import TaskHistoryDB
from divbase_api.worker.worker_db import SyncSessionLocal
from divbase_lib.queries import TaskHistoryResults

logger = logging.getLogger(__name__)


def get_task_history_list(
    user_id: int,
    is_admin: bool = False,
    display_limit: int = 50,
) -> TaskHistoryResults:
    """
    Get a list of the task history from the Flower API.
    """
    api_limit = min(100, display_limit * 5)
    request_url = f"{settings.flower.url}/api/tasks?limit={api_limit}"
    all_tasks = _make_flower_request(request_url)
    filtered_tasks = {}
    for tid, data in all_tasks.items():
        if _check_if_user_has_permission_to_view_task(
            task_id=tid,
            user_id=user_id,
            is_admin=is_admin,
        ):
            filtered_tasks[tid] = data
    return TaskHistoryResults(tasks=filtered_tasks)


def get_task_history_by_id(
    task_id: str,
    user_id: int,
    is_admin: bool = False,
) -> TaskHistoryResults:
    """
    Get the task history from the Flower API for a specific task ID.
    """
    request_url = f"{settings.flower.url}/api/task/info/{task_id}"
    task = _make_flower_request(request_url)
    if not _check_if_user_has_permission_to_view_task(
        task_id=task_id,
        user_id=user_id,
        is_admin=is_admin,
    ):
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


def _check_if_user_has_permission_to_view_task(task_id: str, user_id: int, is_admin: bool) -> bool:
    if is_admin:
        return True

    with SyncSessionLocal() as session:
        entry = session.query(TaskHistoryDB).filter_by(task_id=task_id).first()
        if entry:
            return entry.user_id == user_id
