import os
from time import sleep

import pytest
import requests
from celery import current_app
from celery.backends.redis import RedisBackend
from kombu.connection import Connection

from divbase_tools.tasks import app, simulate_long_task, simulate_quick_task

FLOWER_URL_TESTING_STACK = (
    "http://localhost:5556"  # TODO: could override this as an env var in the testing compose file
)
flower_user = os.environ.get("FLOWER_USER")
flower_password = os.environ.get("FLOWER_PASSWORD")


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


@pytest.mark.integration
@pytest.mark.parametrize(
    "tasks_to_test, wait_time, expected_queue",
    [
        (simulate_quick_task, 1, "quick"),  # Quick task with a short wait time
        (simulate_long_task, 5, "long"),
    ],
)
def test_task_routing(wait_for_celery_task_completion, tasks_to_test, wait_time, expected_queue):
    """
    This test checks that the task routing is set up correctly for the 'simulate_quick_task' task.
    A worker can be assigned to multiple queues, so the test should check that the task is routed to the correct queue,
    in this case the 'quick' queue. Getting logs on which worker executeted the task is easiest done via the Flower API.

    The logic of the test is as follows:
    1. Get the current task routes and the active queues of the workers.
    2. Create a mapping of queues to workers.
    3. Submit a task that is routed to go to the to the 'quick' queue (defined in tasks.py, registered in the app.conf.task_routes).
    5. Use the Flower API to get the task result and check which worker executed the task.
    6. Assert that the worker is in the list of workers assigned to the 'quick' queue.
    """
    current_task_routes = current_app.conf.task_routes
    queue_info = current_app.control.inspect().active_queues()

    current_queue_to_workers_assignment = {}
    for worker, queues in queue_info.items():
        for queue in queues:
            qname = queue["name"]
            if qname not in current_queue_to_workers_assignment:
                current_queue_to_workers_assignment[qname] = []
            current_queue_to_workers_assignment[qname].append(worker)

    result = tasks_to_test.apply_async(wait_time=wait_time)

    task_id = result.id
    assert current_task_routes["tasks.simulate_quick_task"]

    wait_for_celery_task_completion(task_id=task_id, max_wait=30)
    flower_url = f"{FLOWER_URL_TESTING_STACK}/api/tasks"
    auth = (flower_user, flower_password)

    try:
        response = requests.get(flower_url, auth=auth, timeout=3)
    except requests.exceptions.RequestException as e:
        pytest.fail(f"Could not connect to Flower API: {e}")

    if response.ok:
        worker = response.json().get(task_id, {}).get("worker")
        assert worker in current_queue_to_workers_assignment.get(expected_queue), (
            f"Expected task {task_id} to be executed by a worker in the {expected_queue} queue, but got {worker}"
        )
    else:
        print(f"Task {task_id} not found in Flower.")
        raise AssertionError(f"Task {task_id} not found in Flower. Response: {response.text}")


####
## Helper functions here for now to keep the experiment branch more contained. Should be moved to conftest later.
def get_concurrency_of_worker_containers() -> dict:
    """
    Return the concurrency for each running celery worker.
    In Python, use the `app.control.inspect().stats()` method to get the concurrency of each worker.
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
