import ast
import json
import logging
import pickle
from typing import Any

import httpx
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.api_config import settings
from divbase_api.crud.task_history import (
    check_user_can_view_task_id,
    get_task_ids_for_project,
    get_task_ids_for_user,
    get_task_ids_for_user_and_project,
    get_tasks_by_task_id_pg,
)
from divbase_api.exceptions import AuthorizationError, TaskNotFoundInBackendError
from divbase_lib.api_schemas.queries import BcftoolsQueryKwargs, SampleMetadataQueryKwargs
from divbase_lib.api_schemas.task_history import (
    BcftoolsQueryTaskResult,
    DimensionUpdateTaskResult,
    SampleMetadataQueryTaskResult,
    TaskHistoryResult,
    TaskHistoryResults,
)

logger = logging.getLogger(__name__)


# TODO ideally, the flower API could be queried with a list of tasks, but that seem not to be the case. for now, set a limit of 500 tasks to not put a lot of overhead on the call
API_LIMIT = 500
REQUEST_URL_WITH_LIMIT = f"{settings.flower.url}/api/tasks?limit={API_LIMIT}"

# TODO make a workaround to check if all allowed ids were returned? could call the Flower API task_ID by task_ID... but it would be inefficient


async def get_user_task_history_from_postgres(
    db: AsyncSession,
    user_id: int,
    is_admin: bool = False,
) -> TaskHistoryResults:
    """
    Fetch task history for a user for all their projects.
    """

    allowed_task_ids = await get_task_ids_for_user(db, user_id, is_admin)
    # TODO this could be combined with get_tasks_by_task_id_pg into a single DB lookup
    # TODO this function and get_user_and_project_task_history_postgres are very similar. could consider making it more DRY.
    if not allowed_task_ids:
        return TaskHistoryResults(tasks={})

    celery_tasks = await get_tasks_by_task_id_pg(db, allowed_task_ids)

    filtered_tasks = {}
    for task in celery_tasks:
        deserialized = _deserialize_celery_task_metadata(task)
        filtered_tasks[task["task_id"]] = TaskHistoryResult(**deserialized)

    return TaskHistoryResults(tasks=filtered_tasks)


async def get_user_and_project_task_history_postgres(
    db: AsyncSession,
    user_id: int,
    project_id: int | None = None,
    is_admin: bool = False,
) -> TaskHistoryResults:
    """
    Get a list of the task history from the Flower API for a user and project.

    For the case of a results backend purge (task not in flower API results):
    allowed_task_ids uses a db lookup, but all_tasks is fetched from the Flower API.
    Thus, if a task is purged in the results backend, it is naturally excluded.
    """
    allowed_task_ids = await get_task_ids_for_user_and_project(db, user_id, project_id, is_admin)
    # TODO this could be combined with get_tasks_by_task_id_pg into a single DB lookup
    if not allowed_task_ids:
        return TaskHistoryResults(tasks={})

    celery_tasks = await get_tasks_by_task_id_pg(db, allowed_task_ids)

    filtered_tasks = {}
    for task in celery_tasks:
        deserialized = _deserialize_celery_task_metadata(task)
        filtered_tasks[task["task_id"]] = TaskHistoryResult(**deserialized)

    return TaskHistoryResults(tasks=filtered_tasks)


async def get_project_task_history(
    db: AsyncSession,
    project_id: int,
) -> TaskHistoryResults:
    """
    Get the task history of a project from the Flower API.

    """

    allowed_task_ids = await get_task_ids_for_project(db, project_id)
    # TODO this could be combined with get_tasks_by_task_id_pg into a single DB lookup
    if not allowed_task_ids:
        return TaskHistoryResults(tasks={})

    celery_tasks = await get_tasks_by_task_id_pg(db, allowed_task_ids)

    filtered_tasks = {}
    for task in celery_tasks:
        deserialized = _deserialize_celery_task_metadata(task)
        filtered_tasks[task["task_id"]] = TaskHistoryResult(**deserialized)

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

    user_is_allowed_to_access_task_id = await check_user_can_view_task_id(
        db=db, task_id=task_id, user_id=user_id, is_admin=is_admin
    )

    if not user_is_allowed_to_access_task_id:
        raise AuthorizationError("Task ID not found or you don't have permission to view the history for this task ID.")

    celery_tasks = await get_tasks_by_task_id_pg(db, {task_id})

    if not celery_tasks:
        raise TaskNotFoundInBackendError()

    task = celery_tasks[0]
    deserialized = _deserialize_celery_task_metadata(task)

    return TaskHistoryResults(tasks={task_id: TaskHistoryResult(**deserialized)})


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

    # Workaround: Handle 404 from Flower API (task not found or no tasks exist) as empty dict and let layers above handle this with a custom error
    if response.status_code == 404:
        return {}

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
            task_data = _assign_response_models_to_flower_task_fields(task_data)
            filtered_tasks[tid] = TaskHistoryResult(**task_data)

    return TaskHistoryResults(tasks=filtered_tasks)


def _assign_response_models_to_flower_task_fields(task_data: dict) -> dict:
    """Parse the 'result' field of a Flower task data dictionary into the appropriate model."""
    result_raw = task_data.get("result")
    parsed_result = None
    if result_raw:
        try:
            result_dict = ast.literal_eval(result_raw) if isinstance(result_raw, str) else result_raw
            if task_data.get("name") == "tasks.sample_metadata_query":
                parsed_result = SampleMetadataQueryTaskResult(**result_dict)
            elif task_data.get("name") == "tasks.bcftools_query":
                parsed_result = BcftoolsQueryTaskResult(**result_dict)
            elif task_data.get("name") == "tasks.update_vcf_dimensions_task":
                parsed_result = DimensionUpdateTaskResult(**result_dict)
            else:
                parsed_result = result_dict
        except Exception:
            parsed_result = result_raw

    kwargs_raw = task_data.get("kwargs")
    parsed_kwargs = None
    if kwargs_raw:
        try:
            kwargs_dict = ast.literal_eval(kwargs_raw) if isinstance(kwargs_raw, str) else kwargs_raw
            if task_data.get("name") == "tasks.sample_metadata_query":
                parsed_kwargs = SampleMetadataQueryKwargs(**kwargs_dict)
            elif task_data.get("name") == "tasks.bcftools_query":
                parsed_kwargs = BcftoolsQueryKwargs(**kwargs_dict)
            else:
                parsed_kwargs = kwargs_dict
        except Exception:
            parsed_kwargs = kwargs_raw

    task_data["result"] = parsed_result
    task_data["kwargs"] = parsed_kwargs
    return task_data


def _deserialize_celery_task_metadata(task: dict) -> dict:
    """
    Helper function to deserialize tasks from the db CeleryTaskMeta table (SQLalchemy+postgres celery results backend).
    The results backend serializes fields controlled by result_extended=True (e.g. args, kwargs) as JSON, but the task results
    as pickle.

    Handles task args, although the pattern for divbase celery tasks is to use kwargs.
    """

    args = []
    if task.get("args"):
        try:
            args_bytes = task["args"]
            args_str = args_bytes.decode("utf-8") if isinstance(args_bytes, bytes) else args_bytes
            args = json.loads(args_str)
        except (json.JSONDecodeError, UnicodeDecodeError) as e:
            logger.warning(f"Failed to decode args for task {task.get('task_id')}: {e}")

    kwargs = {}
    if task.get("kwargs"):
        try:
            kwargs_bytes = task["kwargs"]
            kwargs_str = kwargs_bytes.decode("utf-8") if isinstance(kwargs_bytes, bytes) else kwargs_bytes
            kwargs = json.loads(kwargs_str)
        except (json.JSONDecodeError, UnicodeDecodeError) as e:
            logger.warning(f"Failed to decode kwargs for task {task.get('task_id')}: {e}")

    result_data = {}
    if task.get("result"):
        result_bytes = task["result"]
        if isinstance(result_bytes, bytes) and result_bytes[:1] == b"\x80":
            try:
                result_data = pickle.loads(result_bytes)
            except Exception as e:
                logger.warning(f"Failed to unpickle result for task {task.get('task_id')}: {e}")
        else:
            try:
                result_str = result_bytes.decode("utf-8") if isinstance(result_bytes, bytes) else result_bytes
                result_data = json.loads(result_str)
            except (json.JSONDecodeError, UnicodeDecodeError) as e:
                logger.warning(f"Failed to decode JSON result for task {task.get('task_id')}: {e}")

    parsed_result = result_data
    parsed_kwargs = kwargs
    task_name = task.get("name")

    if task_name == "tasks.sample_metadata_query":
        parsed_result = SampleMetadataQueryTaskResult(**result_data) if result_data else None
        parsed_kwargs = SampleMetadataQueryKwargs(**kwargs) if kwargs else None
    elif task_name == "tasks.bcftools_query":
        parsed_result = BcftoolsQueryTaskResult(**result_data) if result_data else None
        parsed_kwargs = BcftoolsQueryKwargs(**kwargs) if kwargs else None
    elif task_name == "tasks.update_vcf_dimensions_task":
        parsed_result = DimensionUpdateTaskResult(**result_data) if result_data else None

    args_str_for_flower = json.dumps(args) if isinstance(args, list) else str(args)

    runtime = None
    started_at = task.get("started_at")
    completed_at = task.get("completed_at")
    if started_at and completed_at:
        runtime = (completed_at - started_at).total_seconds()

    return {
        "uuid": task.get("task_id"),
        "status": task.get("status"),
        "result": parsed_result,
        "date_done": task.get("date_done").isoformat() if task.get("date_done") else None,
        "name": task_name,
        "args": args_str_for_flower,
        "kwargs": parsed_kwargs,
        "worker": task.get("worker"),
        "created_at": task.get("created_at").timestamp() if task.get("created_at") else None,
        "started_at": started_at.timestamp() if started_at else None,
        "completed_at": completed_at.timestamp() if completed_at else None,
        "runtime": runtime,
    }
