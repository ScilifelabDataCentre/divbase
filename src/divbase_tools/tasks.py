import logging
import os
from time import sleep

from celery import Celery

from divbase_tools.queries import pipe_query_command

CELERY_BROKER_URL = "pyamqp://guest@rabbitmq//"

broker_url = os.environ.get("CELERY_BROKER_URL", "pyamqp://guest@localhost//")
app = Celery("divbase_tools", broker=broker_url, backend="rpc://")

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
    """
    task_id = bcftools_pipe_task.request.id
    logger.info(f"Starting bcftools pipe task with Celery task ID: {task_id}")

    pipe_query_command(command=command, bcftools_inputs=bcftools_inputs)
