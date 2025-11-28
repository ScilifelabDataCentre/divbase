import json
import logging
import pickle

from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.task_history import (
    get_task_by_id_if_user_allowed,
    get_tasks_for_project_pg,
    get_tasks_for_user_and_project_pg,
    get_tasks_for_user_pg,
)
from divbase_api.exceptions import AuthorizationError
from divbase_lib.api_schemas.queries import BcftoolsQueryKwargs, SampleMetadataQueryKwargs
from divbase_lib.api_schemas.task_history import (
    BcftoolsQueryTaskResult,
    DimensionUpdateTaskResult,
    SampleMetadataQueryTaskResult,
    TaskHistoryResult,
    TaskHistoryResults,
)

logger = logging.getLogger(__name__)


# TODO make a workaround to check if all allowed ids were returned? could call the Flower API task_ID by task_ID... but it would be inefficient


async def get_user_task_history_from_postgres(
    db: AsyncSession,
    user_id: int,
    is_admin: bool = False,
) -> TaskHistoryResults:
    """
    Fetch task history for a user for all their projects.
    """

    celery_tasks = await get_tasks_for_user_pg(db, user_id, is_admin)

    if not celery_tasks:
        return TaskHistoryResults(tasks={})

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
    Fetch task history for a user for a specific project they belong to.
    """

    celery_tasks = await get_tasks_for_user_and_project_pg(db, user_id, project_id, is_admin)

    if not celery_tasks:
        return TaskHistoryResults(tasks={})

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

    celery_tasks = await get_tasks_for_project_pg(db, project_id)

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

    celery_task = await get_task_by_id_if_user_allowed(db=db, task_id=task_id, user_id=user_id, is_admin=is_admin)

    if not celery_task:
        raise AuthorizationError("Task ID not found or you don't have permission to view the history for this task ID.")

    deserialized = _deserialize_celery_task_metadata(celery_task)
    return TaskHistoryResults(tasks={task_id: TaskHistoryResult(**deserialized)})


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
