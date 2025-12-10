import json
import logging
import pickle

from divbase_lib.api_schemas.queries import (
    BcftoolsQueryKwargs,
    BcftoolsQueryTaskResult,
    SampleMetadataQueryKwargs,
    SampleMetadataQueryTaskResult,
)
from divbase_lib.api_schemas.task_history import (
    TaskHistoryResult,
    TaskHistoryResults,
)
from divbase_lib.api_schemas.vcf_dimensions import DimensionUpdateKwargs, DimensionUpdateTaskResult

logger = logging.getLogger(__name__)


def deserialize_tasks_to_result(serialized_tasks: list[dict]) -> TaskHistoryResults:
    """
    Convert a list of task dicts from the DB into a TaskHistoryResults object.
    """
    if not serialized_tasks:
        return TaskHistoryResults(tasks={})

    deserialized_tasks = {}
    for task in serialized_tasks:
        deserialized_task = _deserialize_celery_task_metadata(task)
        deserialized_tasks[task["task_id"]] = TaskHistoryResult(**deserialized_task)

    return TaskHistoryResults(tasks=deserialized_tasks)


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

    task_name = task.get("name")

    is_error_result = isinstance(result_data, dict) and (
        "exc_type" in result_data or "exc_message" in result_data or result_data.get("status") == "error"
    )
    if is_error_result:
        parsed_result = result_data
        parsed_kwargs = kwargs
    else:
        if task_name == "tasks.sample_metadata_query":
            parsed_result = SampleMetadataQueryTaskResult(**result_data) if result_data else None
            parsed_kwargs = SampleMetadataQueryKwargs(**kwargs) if kwargs else None
        elif task_name == "tasks.bcftools_query":
            parsed_result = BcftoolsQueryTaskResult(**result_data) if result_data else None
            parsed_kwargs = BcftoolsQueryKwargs(**kwargs) if kwargs else None
        elif task_name == "tasks.update_vcf_dimensions_task":
            parsed_result = DimensionUpdateTaskResult(**result_data) if result_data else None
            parsed_kwargs = DimensionUpdateKwargs(**kwargs) if kwargs else None
        else:
            # Fallback for Unknown task type - keep everything as dicts
            parsed_result = result_data
            parsed_kwargs = kwargs

    args_as_str = json.dumps(args) if isinstance(args, list) else str(args)

    runtime = None
    started_at = task.get("started_at")
    completed_at = task.get("completed_at")
    if started_at and completed_at:
        runtime = (completed_at - started_at).total_seconds()

    return {
        "uuid": task.get("task_id"),
        "submitter_email": task.get("submitter_email"),
        "status": task.get("status"),
        "result": parsed_result,
        "date_done": task.get("date_done").isoformat() if task.get("date_done") else None,
        "name": task_name,
        "args": args_as_str,
        "kwargs": parsed_kwargs,
        "worker": task.get("worker"),
        "created_at": task.get("created_at").timestamp() if task.get("created_at") else None,
        "started_at": started_at.timestamp() if started_at else None,
        "completed_at": completed_at.timestamp() if completed_at else None,
        "runtime": runtime,
    }
