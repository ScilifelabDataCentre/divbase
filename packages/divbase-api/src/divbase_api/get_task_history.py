import logging
from typing import Any

import httpx

from divbase_api.api_config import settings

logger = logging.getLogger(__name__)


def get_task_history(task_id: str = None, display_limit: int = 10) -> list:
    """
    Get the task history from the Flower API.

    TODO - consider split this into two functions - one for single task, one for list of tasks.
    TODO - consider a more defined return type than list, e.g. dataclass called TaskItem with fields for id, state, timestamp.
        This would then involve changes to the CLI code too, which currently parses this nested dict.

    Returns a list of tasks.
    """
    if task_id:
        request_url = f"{settings.flower.url}/api/task/info/{task_id}"
    else:
        # TODO - if multiple users, this will get tasks from all users and not give correct number of results back.
        api_limit = min(
            100, display_limit * 5
        )  # return more tasks than requested by the --limit arg to allow for downstream sorting
        request_url = f"{settings.flower.url}/api/tasks?limit={api_limit}"

    tasks = _make_flower_request(request_url)

    task_items = []
    if task_id:
        task_items = [(task_id, tasks, parse_timestamp(tasks.get("started", 0)))]
    else:
        task_items = [(tid, data, parse_timestamp(data.get("started", 0))) for tid, data in tasks.items()]

    return task_items


def parse_timestamp(timestamp):
    """
    Convert timestamp to float, if it is a numeric string.
    """
    if isinstance(timestamp, str) and timestamp.replace(".", "").isdigit():
        return float(timestamp)
    return timestamp


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
