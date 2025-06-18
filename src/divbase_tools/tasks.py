import logging
import os
from time import sleep

from celery import Celery

from divbase_tools.queries import pipe_query_command

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

logger = logging.getLogger(__name__)


@app.task(name="tasks.add")
def add(x, y):
    """
    A Hello World-like task to learn how to run tasks in Celery.
    """
    i = 0
    while i < 5:
        sleep(1)
        print("processing...")
        i += 1
    return x + y


@app.task(name="tasks.bcftools_pipe")
def bcftools_pipe_task(command, bcftools_inputs):
    """
    Run pipe_query_command as a Celery task.

    The return messages defined here can be fetched by:
    result = app.AsyncResult(task_id); result.get()
    """
    task_id = bcftools_pipe_task.request.id
    logger.info(f"Starting bcftools pipe task with Celery task ID: {task_id}")

    try:
        pipe_query_command(command=command, bcftools_inputs=bcftools_inputs)
        return {"status": "completed", "output_file": "merged.vcf.gz"}
    except Exception as e:
        logger.error(f"Error in bcftools task: {str(e)}")
        return {"status": "error", "error": str(e)}
