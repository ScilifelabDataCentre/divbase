from unittest.mock import MagicMock, patch

from divbase_tools.tasks import bcftools_pipe_task


@patch("divbase_tools.tasks.BcftoolsQueryManager")
def test_bcftools_pipe_task_directly(mock_bcftools_manager, example_sidecar_metadata_inputs_outputs):
    """Run bcftools_pipe_task direcly without celery. Use mocked dependencies."""

    mock_manager_instance = MagicMock()
    mock_bcftools_manager.return_value = mock_manager_instance
    output_file = "merged.vcf.gz"
    mock_manager_instance.execute_pipe.return_value = output_file

    command = "view -s SAMPLES"
    result = bcftools_pipe_task(command, example_sidecar_metadata_inputs_outputs, submitter=None)

    mock_bcftools_manager.assert_called_once()
    mock_manager_instance.execute_pipe.assert_called_once_with(command, example_sidecar_metadata_inputs_outputs)

    assert result == {"status": "completed", "output_file": f"{output_file}", "submitter": None}
