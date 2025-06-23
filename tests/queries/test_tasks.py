from unittest.mock import MagicMock, patch

from divbase_tools.tasks import app, bcftools_pipe_task


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
    result = bcftools_pipe_task(command, example_sidecar_metadata_inputs_outputs, submitter=None)

    mock_bcftools_manager.assert_called_once()
    mock_manager_instance.execute_pipe.assert_called_once_with(command, example_sidecar_metadata_inputs_outputs)

    assert result == {"status": "completed", "output_file": f"{output_file}", "submitter": None}


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

            result = bcftools_pipe_task.delay(command, {"filenames": ["test.vcf"]}, submitter=None)

            task_result = result.get()

            assert task_result["status"] == "completed"
            assert task_result["output_file"] == output_file
    finally:
        app.conf.task_always_eager = original_task_always_eager_value
        app.conf.task_eager_propagates = original_task_eager_propagates_value
