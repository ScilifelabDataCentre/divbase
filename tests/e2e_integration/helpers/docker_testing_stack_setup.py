"""
Helper module to set up and tear down the docker compose stack used for testing purposes.

When running tests in GH actions we don't need to build the images, so we supply an additional override compose file to do this
And update the compose commands accordingly (e.g. drop the --build flag).

When running with --coverage-docker, a third compose overlay instruments FastAPI and
the Celery workers with coverage.py. Containers are stopped gracefully (SIGTERM) so the
coverage atexit handler flushes .coverage.* files before the container exits.
"""

import json
import os
import subprocess
import time
from pathlib import Path

import httpx

DOCKER_COMPOSE_FILE = "docker/divbase_compose.yaml"
DOCKER_COMPOSE_OVERIDE_FILE = "docker/divbase_compose.tests.yaml"
DOCKER_COMPOSE_CI_OVERRIDE = "docker/divbase_compose.tests.ci.yaml"
DOCKER_COMPOSE_COVERAGE_OVERRIDE = "docker/divbase_compose.tests.coverage.yaml"
COVERAGE_DATA_DIR = Path("docker/coverage-data")
TESTING_STACK_NAME = "divbase-tests"

GH_ACTION_RUN = os.getenv("GITHUB_ACTIONS_RUNNER") == "true"


def _compose_base(coverage_mode: bool = False) -> list[str]:
    """Build the docker compose command prefix as a list for subprocess."""
    parts = ["docker", "compose", "-f", DOCKER_COMPOSE_FILE, "-f", DOCKER_COMPOSE_OVERIDE_FILE]
    if GH_ACTION_RUN:
        parts += ["-f", DOCKER_COMPOSE_CI_OVERRIDE]
    if coverage_mode:
        parts += ["-f", DOCKER_COMPOSE_COVERAGE_OVERRIDE]
    return parts


def start_compose_stack(coverage_mode: bool = False) -> None:
    """Start job system docker stack, then wait for all services to become healthy."""
    if coverage_mode:
        _pre_create_coverage_data_dir()

    up_args = ["up", "-d"] if GH_ACTION_RUN else ["up", "-d", "--build"]
    subprocess.run(_compose_base(coverage_mode) + up_args, check=True)
    print("Waiting for job system docker stack to start and perform healthchecks...")
    time.sleep(2)
    wait_for_docker_stack_healthy(stack_name=TESTING_STACK_NAME, timeout=120)
    _wait_for_api_ready()


def stop_compose_stack(coverage_mode: bool = False) -> None:
    """
    Stop job system docker stack.

    In e2e coverage mode (pytest tests/e2e_integration/ --coverage-docker), send SIGTERM first (docker compose stop) so coverage.py atexit
    handlers flush .coverage.* files to docker/coverage-data/ before containers exit. Using -v in the down command ensures anonymous docker
    volumes created by tmpfs are removed from disk upon tear down, preventing disk space leaks.
    """
    base = _compose_base(coverage_mode)

    if coverage_mode:
        subprocess.run(base + ["stop", "--timeout", "30"], check=True)
        print("\nDocker stack stopped gracefully — coverage data flushed to docker/coverage-data/")

    subprocess.run(base + ["down", "-v"], check=True)
    print("\nStopped job system docker stack.")


def _wait_for_api_ready(
    url: str = "http://localhost:8001/api/v1/core/health", retries: int = 20, delay: int = 5
) -> None:
    """Confirm the FastAPI server is actually accepting connections.

    Docker health checks confirm the container is responding internally, but there can be a brief
    window where the external port mapping is not yet ready or the app is still initialising.
    This provides a second layer of assurance before setup_test_data() is called.
    """
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            response = httpx.get(url, timeout=5)
            if response.status_code == 200:
                print(f"API ready at {url}")
                return
            last_error = RuntimeError(f"Unexpected status {response.status_code}")
        except httpx.ConnectError as e:
            last_error = e
        print(f"API not ready yet (attempt {attempt}/{retries}): {last_error}")
        time.sleep(delay)
    raise RuntimeError(f"API at {url} failed to respond after {retries} attempts. Last error: {last_error}")


def _pre_create_coverage_data_dir() -> None:
    """Create docker/coverage-data/ with open permissions before the stack starts.

    Docker would create the directory as root if it doesn't exist, causing permission
    errors when appuser (UID 1000) inside the containers tries to write coverage files.
    """
    COVERAGE_DATA_DIR.mkdir(parents=True, exist_ok=True)
    COVERAGE_DATA_DIR.chmod(0o777)
    print(f"Coverage data directory ready: {COVERAGE_DATA_DIR}")


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
        # --all includes exited containers so we can detect crashed services immediately.
        result = subprocess.run(
            ["docker", "compose", "-p", stack_name, "ps", "--all", "--format", "json"], capture_output=True, text=True
        )
        container_status = [json.loads(line) for line in result.stdout.strip().splitlines() if line.strip()]

        if not container_status:
            print("No containers found yet, waiting...")
            time.sleep(5)
            continue

        healthy = []
        unhealthy = []
        for c in container_status:
            state = c.get("State")
            health = c.get("Health", "")
            exit_code = c.get("ExitCode", 0)
            service = c.get("Service") or c.get("Name")

            if state == "exited":
                if exit_code == 0:
                    # Init container (e.g. db-migrator, minio-setup) completed successfully — skip.
                    pass
                else:
                    raise RuntimeError(f"Service '{service}' exited unexpectedly with code {exit_code}")
            elif state == "running":
                if not health or health == "healthy":
                    healthy.append(c)
                else:
                    unhealthy.append(c)  # health == "starting" or "unhealthy"
            else:
                unhealthy.append(c)

        print(
            "Healthcheck passed (or not configured in compose file):",
            [c.get("Service") or c.get("Name") for c in healthy],
        )
        print("Healthcheck not yet passed:", [c.get("Service") or c.get("Name") for c in unhealthy])

        if not unhealthy:
            print("All containers healthy!")
            return
        print("Waiting for containers to become healthy...")
        time.sleep(5)
    raise RuntimeError("Some containers failed to become healthy in time.")
