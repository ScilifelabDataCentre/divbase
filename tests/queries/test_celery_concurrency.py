from time import sleep

import pytest
from celery import current_app
from celery.backends.redis import RedisBackend
from kombu.connection import Connection

from divbase_tools.tasks import app, simulate_quick_task


@pytest.mark.integration
def test_concurrency_of_worker_containers_connected_to_default_queue(wait_for_celery_task_completion):
    """
    This test checks that multiple tasks can run concurrently on the worker container. Each Celery worker can run multiple tasks concurrently and
    unless specifically set when initiating the celery app, the concurrency of the worker is based on the number of CPUs available on the host machine.
    If there are multiple worker (containers) running, the concurrency is the sum of the concurrency of each worker.

    Then submit multiple tasks to the worker and check that they are running concurrently. Specifically, submit 1 more task than the total concurrency of the worker containers.
    By timing the sleep cycle of the task itself and the how long we wait for the tasks to complete, we can check that multiple tasks are running concurrently and that the one task
    that exceeds the concurrency limit is PENDING until one of the other tasks finishes.

    Note! This test specifically looks for the 'celery' queue, which is the default queue for Celery tasks. Multiple workers can be assigned to the 'celery' queue,
    and the test should be able to accomodate for that case. However, celery does not automatically adjust the max possible concurrency across containers so there is a risk
    that when multiple workers are assigned to the 'celery' queue (without specifying the concurrency parameter in the container), the concurrency will be higher than the CPU count
    of the host machine, which could lead to excessive load on the host machine and potentially cause the test to fail.
    """

    try:
        broker_url = app.conf.broker_url
        with Connection(broker_url) as conn:
            conn.ensure_connection(max_retries=1)

        if isinstance(app.backend, RedisBackend):
            app.backend.client.ping()

    except Exception as e:
        pytest.skip(f"This test requires services not available: {str(e)}. Is Docker Compose running?")

    wait_time = 5
    max_concurrency_per_worker = get_concurrency_of_worker_containers()
    concurrency_of_workers_connected_to_default_queue = get_concurrency_of_default_queue(max_concurrency_per_worker)
    total_concurrency_on_host = sum(concurrency_of_workers_connected_to_default_queue.values())
    task_count = total_concurrency_on_host + 1

    async_results = [simulate_quick_task.apply_async(wait_time=wait_time, queue="celery") for _ in range(task_count)]

    sleep(1)

    states = [result.state for result in async_results]
    print(f"Task states after 1 second: {states}")

    started_count = states.count("STARTED")
    pending_count = states.count("PENDING")
    assert started_count > 1, "Expected multiple tasks to run concurrently"
    assert started_count <= total_concurrency_on_host, (
        f"Expected at most {total_concurrency_on_host} concurrent tasks on the current machine"
    )
    assert pending_count == task_count - started_count, (
        "Expected more tasks to be pending when concurrency limit is exceeded"
    )

    for result in async_results:
        wait_for_celery_task_completion(task_id=result.id, max_wait=30)


def get_concurrency_of_worker_containers() -> dict:
    """
    Return the concurrency for each running celery worker.
    In Python, you can use the `app.control.inspect().stats()` method to get the concurrency of each worker.
    (From the Celery CLI in the terminal, it is: celery -A divbase_tools.tasks inspect stats)
    """

    app_stats = current_app.control.inspect().stats()
    max_concurrency_per_worker = {}
    for worker, stats in app_stats.items():
        max_concurrency_per_worker[worker] = stats.get("pool", {}).get("max-concurrency")

    return max_concurrency_per_worker


def get_concurrency_of_default_queue(max_concurrency_per_worker: dict) -> dict:
    """
    The default queue is the 'celery' queue, which gets used by tasks that do not specify a queue.
    To find which of the workers that are running the 'celery' queue, get the active queues of each worker.
    Then use the worker ID to filter the max_concurrency_per_worker dictionary.
    """

    queue_info = current_app.control.inspect().active_queues()
    for worker, queues in queue_info.items():
        for queue in queues:
            if queue["name"] == "celery":
                return {worker: max_concurrency_per_worker[worker]}
