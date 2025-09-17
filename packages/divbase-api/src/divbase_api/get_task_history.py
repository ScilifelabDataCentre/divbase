import logging
import os

import httpx

logger = logging.getLogger(__name__)

FLOWER_URL = os.environ.get("FLOWER_HOST", "http://localhost:5555")


def get_task_history(task_id: str = None, display_limit: int = 10) -> list:
    """
    Get the task history from the Flower API.

    Returns a list of tasks.
    """

    flower_user = os.environ.get("FLOWER_USER")
    flower_password = os.environ.get("FLOWER_PASSWORD")
    divbase_user = os.environ.get("DIVBASE_USER")

    # TODO - these should be considered required, not warnings.
    if not flower_user:
        logger.warning("FLOWER_USER not provided or set in environment")
    if not flower_password:
        logger.warning("FLOWER_PASSWORD not provided or set in environment")
    if not divbase_user:
        logger.warning("DIVBASE_USER not set in environment")

    if task_id:
        request_url = f"{FLOWER_URL}/api/task/info/{task_id}"
    else:
        # TODO - if multiple users, this will get tasks from all users and not give correct number of results back.
        api_limit = min(
            100, display_limit * 5
        )  # return more tasks than requested by the --limit arg to allow for downstream sorting
        request_url = f"{FLOWER_URL}/api/tasks?limit={api_limit}"

    auth = (flower_user, flower_password)
    with httpx.Client() as client:
        response = client.get(request_url, auth=auth, timeout=3.0)

    if response.status_code != 200:
        print(f"Failed to fetch tasks for task ID {task_id}. Status code: {response.status_code}")

    tasks = response.json()
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
