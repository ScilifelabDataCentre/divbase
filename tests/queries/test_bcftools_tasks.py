from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from celery import current_app
from celery.backends.redis import RedisBackend
from kombu.connection import Connection

from divbase_tools.tasks import bcftools_pipe_task

FIXTURES_DIR = Path(__file__).parent.parent / "fixtures"


@pytest.mark.unit
@patch("divbase_tools.tasks.BcftoolsQueryManager")
@patch("divbase_tools.s3_client.S3FileManager.download_files")
@patch("divbase_tools.s3_client.S3FileManager.upload_files")
def test_bcftools_pipe_task_directly(
    mock_upload_files,
    mock_download_files,
    mock_bcftools_manager,
    bcftools_pipe_kwargs_fixture,
    tmp_path,
    sample_tsv_file,
):
    """
    Run bcftools_pipe_task direcly without going via a Celery worker or any other services in the job system.
    Needs many mocked services to run, so it mainly tests the task logic and flow.

    bcftools_pipe_kwargs_fixture was designed to use with celery, e.g. bcftools_pipe_task.apply_async(kwargs=task_kwargs),
    but the fixture (which is a dict with the task keys) can be used here if unpacked first.
    """

    dummy_vcf = tmp_path / "test.vcf.gz"
    mock_download_files.side_effect = [
        [sample_tsv_file],
        [dummy_vcf],
    ]
    mock_upload_files.return_value = None

    output_file = "merged.vcf.gz"
    mock_manager_instance = MagicMock()
    mock_manager_instance.execute_pipe.return_value = output_file
    mock_bcftools_manager.return_value = mock_manager_instance

    result = bcftools_pipe_task(**bcftools_pipe_kwargs_fixture)

    mock_bcftools_manager.assert_called_once()
    mock_manager_instance.execute_pipe.assert_called_once()
    mock_download_files.assert_called()
    mock_upload_files.assert_called_once()

    user_name = bcftools_pipe_kwargs_fixture.get("user_name")
    assert result == {
        "status": "completed",
        "output_file": f"{output_file}",
        "submitter": user_name,
    }


@pytest.mark.unit
@patch("divbase_tools.tasks.BcftoolsQueryManager")
@patch("divbase_tools.s3_client.S3FileManager.download_files")
@patch("divbase_tools.queries.run_sidecar_metadata_query")
@patch("divbase_tools.s3_client.S3FileManager.upload_files")
def test_bcftools_pipe_task_using_eager_mode(
    mock_upload_files,
    mock_sidecar_query,
    mock_download_files,
    mock_bcftools_manager,
    bcftools_pipe_kwargs_fixture,
    sample_tsv_file,
    tmp_path,
):
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

    original_task_always_eager_value = current_app.conf.task_always_eager
    original_task_eager_propagates_value = current_app.conf.task_eager_propagates

    try:
        current_app.conf.update(
            task_always_eager=True,
            task_eager_propagates=True,
        )

        dummy_vcf = tmp_path / "test.vcf.gz"
        mock_download_files.side_effect = [
            [sample_tsv_file],
            [dummy_vcf],
        ]
        mock_upload_files.return_value = None

        output_file = "merged.vcf.gz"
        mock_manager_instance = MagicMock()
        mock_manager_instance.execute_pipe.return_value = output_file
        mock_bcftools_manager.return_value = mock_manager_instance

        mock_metadata_result = MagicMock()
        mock_metadata_result.unique_filenames = [str(dummy_vcf)]
        mock_metadata_result.sample_and_filename_subset = [{"SampleID": "S1", "Filename": str(dummy_vcf)}]
        mock_metadata_result.unique_sample_ids = ["S1"]
        mock_sidecar_query.return_value = mock_metadata_result

        result = bcftools_pipe_task.delay(**bcftools_pipe_kwargs_fixture)

        task_result = result.get()

        user_name = bcftools_pipe_kwargs_fixture.get("user_name")
        assert task_result == {
            "status": "completed",
            "output_file": f"{output_file}",
            "submitter": user_name,
        }
    finally:
        current_app.conf.task_always_eager = original_task_always_eager_value
        current_app.conf.task_eager_propagates = original_task_eager_propagates_value


@pytest.mark.integration
def test_bcftools_pipe_task_with_real_worker(wait_for_celery_task_completion, bcftools_pipe_kwargs_fixture):
    """
    Integration test in which bcftools_pipe_task is run with a real Celery worker.
    Runs locally using the docker compose testing stack defined and performs a real sidecar and bcftools
    query by loading VCF files from the tests/fixtures dir. It was designed for having RabbitMQ as the broker, Redis as the backend,
    and a custom Celery worker image that has bcftools installed.
    (this test does not download fixture from bucket, since that is handled by the CLI layer)

    Does not assert file download and uploads, since that is handled by tests in tests/cli_commands/test_query_cli.py.
    (but the current version of the task does I/O with files from the bucket, so it is tested indirectly)

    This test is marked with 'integration' so it can be skipped by: pytest -m "not integration"
    Likewise, all tests marked thusly can be run with: pytest -m integration
    """

    try:
        broker_url = current_app.conf.broker_url
        with Connection(broker_url) as conn:
            conn.ensure_connection(max_retries=1)

        if isinstance(current_app.backend, RedisBackend):
            current_app.backend.client.ping()

    except Exception as e:
        pytest.skip(f"This test requires services not available: {str(e)}. Is Docker Compose running?")

    async_result = bcftools_pipe_task.apply_async(kwargs=bcftools_pipe_kwargs_fixture)
    task_id = async_result.id
    task_result = wait_for_celery_task_completion(task_id=task_id, max_wait=30)

    user_name = bcftools_pipe_kwargs_fixture.get("user_name")
    assert task_result["status"] == "completed"
    assert task_result["submitter"] == user_name
    assert task_result["output_file"].startswith("merged_") and task_result["output_file"].endswith(".vcf.gz")
