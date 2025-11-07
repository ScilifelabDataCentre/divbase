import ast
import logging
from typing import Any

import httpx

from divbase_api.api_config import settings
from divbase_lib.queries import TaskHistoryResults

logger = logging.getLogger(__name__)


def get_task_history_list(
    display_limit: int = 50,
    submitter_email: str = None,
    is_admin: bool = False,
) -> TaskHistoryResults:
    """
    Get a list of the task history from the Flower API.
    """
    api_limit = min(100, display_limit * 5)
    request_url = f"{settings.flower.url}/api/tasks?limit={api_limit}"
    all_tasks = _make_flower_request(request_url)
    filtered_tasks = {
        tid: data
        for tid, data in all_tasks.items()
        if _check_if_user_has_permission_to_view_task(data, submitter_email, is_admin)
    }
    return TaskHistoryResults(tasks=filtered_tasks)


def get_task_history_by_id(
    task_id: str = None,
    submitter_email: str = None,
    is_admin: bool = False,
) -> TaskHistoryResults:
    """
    Get the task history from the Flower API for a specific task ID.
    """
    request_url = f"{settings.flower.url}/api/task/info/{task_id}"
    task = _make_flower_request(request_url)
    if not _check_if_user_has_permission_to_view_task(task, submitter_email, is_admin):
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


def _check_if_user_has_permission_to_view_task(task: dict, submitter_email: str, is_admin: bool) -> bool:
    if is_admin or not submitter_email:
        return True
    kwargs = task.get("kwargs", "{}")
    try:
        parsed_kwargs = ast.literal_eval(kwargs)
    except Exception:
        parsed_kwargs = {}
    submitter = parsed_kwargs.get("user_name")
    return submitter == submitter_email
