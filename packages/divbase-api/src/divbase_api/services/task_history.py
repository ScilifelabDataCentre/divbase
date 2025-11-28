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
from divbase_api.models.task_history import CeleryTaskMeta
from divbase_lib.api_schemas.queries import BcftoolsQueryKwargs, SampleMetadataQueryKwargs
from divbase_lib.api_schemas.task_history import (
    BcftoolsQueryTaskResult,
    DimensionUpdateTaskResult,
    FlowerTaskResult,
    SampleMetadataQueryTaskResult,
    TaskHistoryResults,
)

logger = logging.getLogger(__name__)


# TODO ideally, the flower API could be queried with a list of tasks, but that seem not to be the case. for now, set a limit of 500 tasks to not put a lot of overhead on the call
API_LIMIT = 500
REQUEST_URL_WITH_LIMIT = f"{settings.flower.url}/api/tasks?limit={API_LIMIT}"

# TODO make a workaround to check if all allowed ids were returned? could call the Flower API task_ID by task_ID... but it would be inefficient


async def get_user_task_history(
    db: AsyncSession,
    user_id: int,
    is_admin: bool = False,
) -> TaskHistoryResults:
    """
    Get a list of the task history from the Flower API.

    For the case of a results backend purge (task not in flower API results):
    allowed_task_ids uses a db lookup, but all_tasks is fetched from the Flower API.
    Thus, if a task is purged in the results backend, it is naturally excluded.
    """

    allowed_task_ids = await get_task_ids_for_user(db, user_id, is_admin)

    if not allowed_task_ids:
        return TaskHistoryResults(tasks={})

    all_tasks = _make_flower_request(REQUEST_URL_WITH_LIMIT)

    filtered_results = _filter_flower_results_by_allowed_task_ids(
        all_tasks=all_tasks, allowed_task_ids=allowed_task_ids
    )
    return filtered_results


async def get_user_task_history_from_postgres(
    db: AsyncSession,
    user_id: int,
    is_admin: bool = False,
) -> TaskHistoryResults:
    """Version of get_user_task_history that queries the celery pg backend table"""

    allowed_task_ids = await get_task_ids_for_user(db, user_id, is_admin)

    if not allowed_task_ids:
        return TaskHistoryResults(tasks={})

    celery_tasks = await get_tasks_by_task_id_pg(db, allowed_task_ids)

    filtered_tasks = {}
    for task in celery_tasks:
        deserialized = _deserialize_celery_task_metadata(task)
        filtered_tasks[task.task_id] = FlowerTaskResult(**deserialized)

    return TaskHistoryResults(tasks=filtered_tasks)


async def get_user_and_project_task_history(
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
    Get the task history of a project from the Flower API.

    """

    allowed_task_ids = await get_task_ids_for_project(db, project_id)

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

    user_is_allowed_to_access_task_id = await check_user_can_view_task_id(
        db=db, task_id=task_id, user_id=user_id, is_admin=is_admin
    )

    if not user_is_allowed_to_access_task_id:
        raise AuthorizationError("Task ID not found or you don't have permission to view the history for this task ID.")

    request_url = f"{settings.flower.url}/api/task/info/{task_id}"
    task_data = _make_flower_request(request_url)

    if not task_data:
        raise TaskNotFoundInBackendError()

    task_data = _assign_response_models_to_flower_task_fields(task_data)
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
            filtered_tasks[tid] = FlowerTaskResult(**task_data)

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


def _deserialize_celery_task_metadata(task: CeleryTaskMeta) -> dict:
    """
    Helper function to deserialize tasks from the db CeleryTaskMeta table (SQLalchemy+postgres celery results backend).
    The results backend serializes fields controlled by result_extended=True (e.g. args, kwargs) as JSON, but the task results
    as pickle.

    Handles task args, although the pattern for divbase celery tasks is to use kwargs.
    """

    args = []
    if task.args:
        try:
            args_str = task.args.decode("utf-8") if isinstance(task.args, bytes) else task.args
            args = json.loads(args_str)
        except (json.JSONDecodeError, UnicodeDecodeError) as e:
            logger.warning(f"Failed to decode args for task {task.task_id}: {e}")

    kwargs = {}
    if task.kwargs:
        try:
            kwargs_str = task.kwargs.decode("utf-8") if isinstance(task.kwargs, bytes) else task.kwargs
            kwargs = json.loads(kwargs_str)
        except (json.JSONDecodeError, UnicodeDecodeError) as e:
            logger.warning(f"Failed to decode kwargs for task {task.task_id}: {e}")

    result_data = {}
    if task.result:
        if isinstance(task.result, bytes) and task.result[:1] == b"\x80":
            try:
                result_data = pickle.loads(task.result)
            except Exception as e:
                logger.warning(f"Failed to unpickle result for task {task.task_id}: {e}")
        else:
            try:
                result_str = task.result.decode("utf-8") if isinstance(task.result, bytes) else task.result
                result_data = json.loads(result_str)
            except (json.JSONDecodeError, UnicodeDecodeError) as e:
                logger.warning(f"Failed to decode JSON result for task {task.task_id}: {e}")

    parsed_result = result_data
    parsed_kwargs = kwargs
    args_str_for_flower = json.dumps(args) if isinstance(args, list) else str(args)

    if task.name == "tasks.sample_metadata_query":
        parsed_result = SampleMetadataQueryTaskResult(**result_data) if result_data else None
        parsed_kwargs = SampleMetadataQueryKwargs(**kwargs) if kwargs else None
    elif task.name == "tasks.bcftools_query":
        parsed_result = BcftoolsQueryTaskResult(**result_data) if result_data else None
        parsed_kwargs = BcftoolsQueryKwargs(**kwargs) if kwargs else None
    elif task.name == "tasks.update_vcf_dimensions_task":
        parsed_result = DimensionUpdateTaskResult(**result_data) if result_data else None

    return {
        "uuid": task.task_id,
        "status": task.status,
        "result": parsed_result,
        "date_done": task.date_done.isoformat() if task.date_done else None,
        "name": task.name,
        "args": args_str_for_flower,
        "kwargs": parsed_kwargs,
        "worker": task.worker,
    }
