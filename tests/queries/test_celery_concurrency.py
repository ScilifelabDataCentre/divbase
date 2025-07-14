import os
from time import sleep

import pytest
import requests
from celery import current_app
from celery.backends.redis import RedisBackend
from kombu.connection import Connection

from divbase_tools.tasks import app, bcftools_pipe_task, sample_metadata_query_task

FLOWER_URL_TESTING_STACK = (
    "http://localhost:5556"  # TODO: could override this as an env var in the testing compose file
)


@pytest.mark.integration
def test_concurrency_of_worker_containers_connected_to_default_queue(
    wait_for_celery_task_completion,
    concurrency_of_default_queue,
    bcftools_pipe_kwargs_fixture,
):
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

    NOTE! This test is potentially flaky since it is dependent on timings. If the completion time for the task is less than 0.01 (which does not seem to be the case during the local dev so far), this test will likely fail.
    """

    try:
        broker_url = app.conf.broker_url
        with Connection(broker_url) as conn:
            conn.ensure_connection(max_retries=1)

        if isinstance(app.backend, RedisBackend):
            app.backend.client.ping()

    except Exception as e:
        pytest.skip(f"This test requires services not available: {str(e)}. Is Docker Compose running?")

    total_default_queue_concurrency_on_host = sum(concurrency_of_default_queue.values())
    task_count = total_default_queue_concurrency_on_host + 1

    async_results = [
        bcftools_pipe_task.apply_async(kwargs=bcftools_pipe_kwargs_fixture, queue="celery") for _ in range(task_count)
    ]

    wait_time = 0.01
    sleep(wait_time)

    states = [result.state for result in async_results]

    started_count = states.count("STARTED")
    pending_count = states.count("PENDING")
    assert started_count > 1, "Expected multiple tasks to run concurrently"
    assert started_count <= total_default_queue_concurrency_on_host, (
        f"Expected at most {total_default_queue_concurrency_on_host} concurrent tasks on the current machine"
    )
    assert pending_count == task_count - started_count, (
        "Expected more tasks to be pending when concurrency limit is exceeded"
    )

    for result in async_results:
        wait_for_celery_task_completion(task_id=result.id, max_wait=30)


@pytest.mark.integration
@pytest.mark.parametrize(
    "tasks_to_test, kwargs_fixture, expected_queue",
    [
        (
            sample_metadata_query_task,
            "sample_metadata_query_kwargs_fixture",
            "quick",
        ),
        (
            bcftools_pipe_task,
            "bcftools_pipe_kwargs_fixture",
            "long",
        ),
    ],
)
def test_task_routing(wait_for_celery_task_completion, tasks_to_test, kwargs_fixture, expected_queue, request):
    """
    This test checks that the task routing is set up correctly for the tasks in the test parameter.
    A worker can be assigned to multiple queues, so the test should check that each task is routed to the correct queue (either via static or dynamic task routing).
    Getting logs on which worker executeted the task is easiest done via the Flower API.

    The logic of the test is as follows:
    1. Get the current task routes and assert that the task to be tested by the current test parameter is in the current task routes.
       The logic is different for static and dynamic routing (see also more details below):
       - If static routing is used, the task should be in the current task routes dictionary.
       - If dynamic routing is used, the current task routes is a tuple of functions that each can be called to return a dictionary with routing information.
    2. Get the active queues of the workers and create a mapping of queues to workers.
    3. Submit a task that is routed to go to the to the queue defined in the router in tasks.py (and registered in the app.conf.task_routes).
    4. Use the Flower API to get the task result and check which worker executed the task and assert that the worker is in the list of workers assigned to the
       expected_queue of the current test parameter

    More details on the logic for asserting static vs dynamic routing for tasks in step 1:
    - if static routing is used, current_task_routes is a dict. the functions in the test parameter tasks_to_test are celery tasks and thus have a .name attribute
    - If dynamic routing is used, current_task_routes is a tuple of functions
      There can be multiple dynamic router functions in the tuple, and each return value from each dynamic router function
      becomes an function object that can be called to route a task.
      each function object folow the this pattern from the Celery router docs:
      'route_task(name, args, kwargs, options, task=None, **kw)'
      By calling the function object with the task name and other parameters, a dictionary with the routing information is
      returned (comparable to the dict of the static routing).

    Pytest does not allow fixtures to be used as parameters in parametrized tests. The workaround is call the fixture inside the test using request.getfixturevalue;
    request is a pytest built-in fixture that allows the value of a fixtures to be accessed by calling the fixture name.
    """

    try:
        broker_url = app.conf.broker_url
        with Connection(broker_url) as conn:
            conn.ensure_connection(max_retries=1)

        if isinstance(app.backend, RedisBackend):
            app.backend.client.ping()

    except Exception as e:
        pytest.skip(f"This test requires services not available: {str(e)}. Is Docker Compose running?")

    ## Step 1
    current_task_routes = current_app.conf.task_routes

    if isinstance(current_task_routes, dict):  # For when static routing is used in the celery app
        assert tasks_to_test.name in current_task_routes
    if isinstance(current_task_routes, tuple):  # For when dynamic routing is used in the celery app
        route_dict = None
        for router_function_object in current_task_routes:
            route_dict = router_function_object(
                name=tasks_to_test.name, args=None, kwargs={}, options={}, task=tasks_to_test
            )
            if route_dict:
                break
        assert route_dict is not None, f"No route found for task {tasks_to_test.name} using dynamic router"
        assert route_dict.get("queue") == expected_queue, (
            f"Dynamic router did not route {tasks_to_test.name} to expected queue '{expected_queue}', got '{route_dict.get('queue')}'"
        )

    ## Step 2
    queue_info = current_app.control.inspect().active_queues()

    current_queue_to_workers_assignment = {}
    for worker, queues in queue_info.items():
        for queue in queues:
            qname = queue["name"]
            if qname not in current_queue_to_workers_assignment:
                current_queue_to_workers_assignment[qname] = []
            current_queue_to_workers_assignment[qname].append(worker)

    ## Step 3
    task_kwargs = request.getfixturevalue(kwargs_fixture)
    result = tasks_to_test.apply_async(kwargs=task_kwargs)
    task_id = result.id

    ## Step 4
    wait_for_celery_task_completion(task_id=task_id, max_wait=30)

    flower_user = os.environ.get("FLOWER_USER")
    flower_password = os.environ.get("FLOWER_PASSWORD")
    flower_base_url = os.environ.get("FLOWER_BASE_URL")

    flower_url = f"{flower_base_url}/api/tasks"
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
        raise AssertionError(f"Task {task_id} not found in Flower. Response: {response.text}")
