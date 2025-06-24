import os
import time
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from celery.backends.redis import RedisBackend
from kombu.connection import Connection

from divbase_tools.cli_commands.query_cli import pipe_query
from divbase_tools.tasks import app, bcftools_pipe_task

FIXTURES_DIR = Path(__file__).parent.parent / "fixtures"


@pytest.mark.unit
@patch("divbase_tools.tasks.BcftoolsQueryManager")
def test_bcftools_pipe_task_directly(mock_bcftools_manager, example_sidecar_metadata_inputs_outputs):
    """
    Run bcftools_pipe_task direcly without celery. Use mocked dependencies.
    This tests tests the buissiness logic of the task without needing a Celery worker.
    """

    mock_manager_instance = MagicMock()
    mock_bcftools_manager.return_value = mock_manager_instance
    output_file = "merged.vcf.gz"
    mock_manager_instance.execute_pipe.return_value = output_file

    command = "view -s SAMPLES"
    result = bcftools_pipe_task(command, example_sidecar_metadata_inputs_outputs, submitter="test_user")

    mock_bcftools_manager.assert_called_once()
    mock_manager_instance.execute_pipe.assert_called_once_with(command, example_sidecar_metadata_inputs_outputs)

    assert result == {
        "status": "completed",
        "output_file": f"{output_file}",
        "submitter": "test_user",
    }


@pytest.mark.unit
def test_bcftools_pipe_task_using_eager_mode(example_sidecar_metadata_inputs_outputs):
    """
    Test Celery task in eager mode (no worker needed).
    This is similar to running the task directly, but it uses Celery's eager mode
    which tests the task execution flow as if it were run by a worker.
    Specifically, it tests the .delay() and .get() methods of the task.

    task_always_eager=True means that celery will run the tasks synchronously,
    as if they were called directly, but still using the task infrastructure.
    Use try/finally to restore the original settings after the test even if the test breaks.

    task_eager_propagates=True, means that exceptions raised in the task will propagate
    to the caller, which is useful for troubleshooting.
    """

    original_task_always_eager_value = app.conf.task_always_eager
    original_task_eager_propagates_value = app.conf.task_eager_propagates

    try:
        app.conf.update(
            task_always_eager=True,
            task_eager_propagates=True,
        )

        with patch("divbase_tools.tasks.BcftoolsQueryManager") as mock_bcftools_manager:
            mock_manager_instance = MagicMock()
            mock_bcftools_manager.return_value = mock_manager_instance
            output_file = "merged.vcf.gz"
            mock_manager_instance.execute_pipe.return_value = output_file
            command = "view -s SAMPLES"

            result = bcftools_pipe_task.delay(command, {"filenames": ["test.vcf"]}, submitter="test_user")

            task_result = result.get()

            assert task_result == {
                "status": "completed",
                "output_file": f"{output_file}",
                "submitter": "test_user",
            }
    finally:
        app.conf.task_always_eager = original_task_always_eager_value
        app.conf.task_eager_propagates = original_task_eager_propagates_value


@pytest.mark.integration
def test_bcftools_pipe_task_with_real_worker(demo_sidecar_metadata_inputs_outputs):
    """
    Integration test in which bcftools_pipe_task is run with a real Celery worker.
    Runs locally using the docker-compose setup defined in the tests/queries/docker-compose.yml file.
    It runs a real bcftools query by loading VCF files from the tests/fixtures dir.
    (this test does not download fixture from bucket, since that is handled by the CLI layer)

    This test requires the setup defined in the docker-compose.yml file for queries.
    It was designed for having RabbitMQ as the broker, Redis as the backend,
    and a custom Celery worker image that has bcftools installed.

    At the time of writing, this test imports the broker and backend URLs from divbase_tools.tasks.
    These match the docker-compose setup. For future refernce, these two calls for broker and backend
    could also be set as a pytest fixture that returns:
    broker_url = os.environ.get("CELERY_BROKER_URL", "pyamqp://guest@localhost//")
    result_backend = os.environ.get("CELERY_RESULT_BACKEND", "redis://localhost:6379/0")

    This test is marked with 'integration' so it can be skipped by: pytest -m "not integration"
    Likewise, all tests marked thusly can be run with: pytest -m integration
    """
    command = "view -s SAMPLES; view -r 21:15000000-25000000"
    output_file = "merged.vcf.gz"

    try:
        broker_url = app.conf.broker_url
        with Connection(broker_url) as conn:
            conn.ensure_connection(max_retries=1)

        if isinstance(app.backend, RedisBackend):
            app.backend.client.ping()

    except Exception as e:
        pytest.skip(f"This test requires services not available: {str(e)}. Is Docker Compose running?")

    async_result = bcftools_pipe_task.delay(command, demo_sidecar_metadata_inputs_outputs, submitter="test_user")

    max_wait = 30
    start_time = time.time()

    while not async_result.ready():  # this data is expecte to take <30s to query, so we can allow us to wait for it
        if time.time() - start_time > max_wait:
            pytest.fail(f"Task timed out after {max_wait} seconds")
        time.sleep(1)

    result = async_result.get()

    assert result == {
        "status": "completed",
        "output_file": f"{output_file}",
        "submitter": "test_user",
    }


@pytest.mark.integration
@patch("divbase_tools.cli_commands.query_cli.download_files_command")
def test_pipe_query_e2e_run_async_false(
    mock_download, demo_sidecar_metadata_inputs_outputs, copy_fixtures_to_mock_download_from_bucket
):
    """
    End-to-end test for the pipe_query function, i.e. the CLI command that runs the bcftools query.
    """

    output_file = Path("merged.vcf.gz")
    if output_file.exists():
        print("Output file exists; deleting it to ensure a clean test run")
        output_file.unlink()

    demo_sidecar_metadata_inputs_outputs["filenames"] = [
        f.replace("/app", ".") for f in demo_sidecar_metadata_inputs_outputs["filenames"]
    ]

    mock_download.side_effect = copy_fixtures_to_mock_download_from_bucket(
        demo_sidecar_metadata_inputs_outputs["filenames"]
    )
    os.environ["DIVBASE_USER"] = "test_user"

    tsv_filter = "Area:West of Ireland,Northern Portugal;Sex:F"
    command = "view -s SAMPLES; view -r 21:15000000-25000000"

    pipe_query(
        tsv_file=Path("tests/fixtures/sample_metadata.tsv"),
        tsv_filter=tsv_filter,
        command=command,
        bucket_name="test-bucket",
        config_path=None,
        run_async=False,
    )

    assert output_file.exists(), "Output file was not created"


# TODO: write a e2e test for the CLI entry, with run_async=True. It should use the MinIO test image to handle file downloads instead of the mock function
