import os
from time import sleep

from celery import Celery

CELERY_BROKER_URL = "pyamqp://guest@rabbitmq//"

broker_url = os.environ.get("CELERY_BROKER_URL", "pyamqp://guest@localhost//")
app = Celery("divbase_tools", broker=broker_url, backend="rpc://")


@app.task(name="tasks.add")
def add(x, y):
    i = 0
    while i < 5:
        sleep(1)
        print("processing...")
        i += 1
    return x + y
