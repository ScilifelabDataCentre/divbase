"""
Helper module to set up and tear down a Minio instance used for testing purposes.
"""

import json
import shlex
import subprocess
import time

DOCKER_COMPOSE_FILE = "docker/job_system_compose.yaml"
DOCKER_COMPOSE_OVERIIDE_FILE = "docker/job_system_compose.tests.yaml"
COMPOSE_COMMAND = shlex.split(f"docker compose -f {DOCKER_COMPOSE_FILE} -f {DOCKER_COMPOSE_OVERIIDE_FILE} up -d")

STOP_COMMAND = shlex.split(f"docker compose -f {DOCKER_COMPOSE_FILE} -f {DOCKER_COMPOSE_OVERIIDE_FILE} down")

TESTING_STACK_NAME = "divbase-job-system-tests"

# FLOWER_HEALTH_CHECK_COMMAND = shlex.split("curl -I http://localhost:5556")


def start_job_system() -> None:
    """Start job system docker stack using Docker compose, the call helper function to ensure that all services in stack are healthy."""
    subprocess.run(COMPOSE_COMMAND, check=True)
    print("Waiting for job system docker stack to start and perform healthchecks...")
    time.sleep(2)
    wait_for_docker_stack_healthy(stack_name=TESTING_STACK_NAME, timeout=120)


def stop_job_system() -> None:
    """Stop job system docker stack."""
    subprocess.run(STOP_COMMAND, check=True)
    print("\nStopping job system docker stack...")


def wait_for_docker_stack_healthy(stack_name, timeout=120):
    """
    Wait until all containers in the docker stack are healthy.

    This function checks the status of all containers in the specified docker stack.
    It is intended to be used with services that have health checks defined in their Docker Compose configuration,
    but it handles services without health checks as well (see 2a below).

    Health status logic:
    1.  For a container to be considered healthy, it must be in the "running" state.
    2a. If a container does NOT have a healthcheck defined, the "Health" field is usually None, an empty string, or not present.
        Such containers are considered healthy as long as their "State" is "running".
    2b. If a container HAS a healthcheck defined, the "Health" field will be:
        - "starting" while the healthcheck is initializing or running.
        - "healthy" once the healthcheck has passed.
        - "unhealthy" if the healthcheck fails.
        Only containers with "Health" set to "healthy" are considered healthy in this case.
    """
    start = time.time()
    while time.time() - start < timeout:
        result = subprocess.run(
            ["docker", "compose", "-p", stack_name, "ps", "--format", "json"], capture_output=True, text=True
        )
        container_status = [json.loads(line) for line in result.stdout.strip().splitlines() if line.strip()]

        healthy = []
        unhealthy = []
        for c in container_status:
            state = c.get("State")
            health = c.get("Health")
            if state == "running" and (not health or health == "healthy"):
                healthy.append(c)
            else:
                unhealthy.append(c)

        print("Healthcheck passed:", [c.get("Service") or c.get("Name") for c in healthy])
        print("Healthcheck not yet passed:", [c.get("Service") or c.get("Name") for c in unhealthy])

        if not unhealthy:
            print("All containers healthy!")
            return
        print("Waiting for containers to become healthy...")
        time.sleep(5)
    raise RuntimeError("Some containers failed to become healthy in time.")
