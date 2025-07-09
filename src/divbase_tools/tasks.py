import logging
import os
from pathlib import Path
from time import sleep

from celery import Celery

from divbase_tools.queries import BcftoolsQueryManager
from divbase_tools.s3_client import create_s3_file_manager

logger = logging.getLogger(__name__)

broker_url = os.environ.get("CELERY_BROKER_URL", "pyamqp://guest@localhost//")
result_backend = os.environ.get("CELERY_RESULT_BACKEND", "redis://localhost:6379/0")

app = Celery("divbase_tools", broker=broker_url, backend=result_backend)

# Redis-specific config
app.conf.update(
    result_expires=3600,
    task_track_started=True,
    task_serializer="json",
    accept_content=["json"],
    result_serializer="json",
)


def dynamic_router(name, args, kwargs, options, task=None, **kw):
    """
    Function to dynamically route tasks based on their name and arguments.
    The input arguments of the function are taken from the Celery docs.

    Advanced logic can be implemented here to route tasks to different queues
    based on the task names or arguments.

    For now, just use a simple routing based on substrings in the task name.

    To use static routing instead, use the following the call instead of the dynamic_router function:
        app.conf.task_routes = {
        "tasks.simulate_quick_task": {"queue": "quick"},
        "tasks.simulate_long_task": {"queue": "long"},
        "tasks.bcftools_pipe": {"queue": "long"},
        }
    """
    if "quick" in name:
        return {"queue": "quick"}
    if "long" in name:
        return {"queue": "long"}
    if name == "tasks.bcftools_pipe":
        return {"queue": "long"}
    return {"queue": "celery"}


app.conf.task_routes = (dynamic_router,)


@app.task(name="tasks.simulate_quick_task")
def simulate_quick_task(wait_time: int = 2):
    """
    A simple task to simulate a quick operation.
    Intended for testing of celery concurrency and task management.
    """
    task_id = simulate_quick_task.request.id
    logger.info(f"Starting quick task with Celery task ID: {task_id}")

    sleep(wait_time)
    return {"status": "completed", "message": "Quick task completed successfully."}


@app.task(name="tasks.simulate_long_task")
def simulate_long_task(wait_time: int = 20):
    """
    A simple task to simulate a long(er) operation.
    Intended for testing of celery concurrency and task management.
    Code duplication with simulate_quick_task is intentional to be able test different
    task routing for the two different tasks.
    """
    task_id = simulate_long_task.request.id
    logger.info(f"Starting long task with Celery task ID: {task_id}")

    sleep(wait_time)
    return {"status": "completed", "message": "Long task completed successfully."}


@app.task(name="tasks.bcftools_pipe")
def bcftools_pipe_task(
    command: str,
    bcftools_inputs: dict,
    submitter: str | None,
    bucket_name: str | None = None,
):
    """
    Run pipe_query_command as a Celery task.

    TODO - bucket_name is currently optional so query_cli tests works, but will later be required.

    The return messages defined here can be fetched by:
    result = app.AsyncResult(task_id); result.get()
    """
    task_id = bcftools_pipe_task.request.id
    logger.info(f"Starting bcftools pipe task with Celery task ID: {task_id}")

    # TODO - move dload of metadata file + filtering of tsv here, then download of files from S3 if needed.

    try:
        output_file = BcftoolsQueryManager().execute_pipe(command, bcftools_inputs)
    except Exception as e:
        logger.error(f"Error in bcftools task: {str(e)}")
        return {"status": "error", "error": str(e), "submitter": submitter}

    if bucket_name:
        upload_results_file(output_file=Path(output_file), bucket_name=bucket_name)

    return {"status": "completed", "output_file": f"{output_file}", "submitter": submitter}


def upload_results_file(output_file: Path, bucket_name: str) -> None:
    """
    Upon completion of the task, upload the results file to the specified bucket.
    """
    s3_file_manager = create_s3_file_manager()
    _ = s3_file_manager.upload_files(
        to_upload={output_file.name: output_file},
        bucket_name=bucket_name,
    )
