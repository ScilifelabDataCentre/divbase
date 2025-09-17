import subprocess
from unittest.mock import MagicMock, patch

import pytest

from divbase_lib.exceptions import (
    BcftoolsCommandError,
    BcftoolsEnvironmentError,
    BcftoolsPipeEmptyCommandError,
    BcftoolsPipeUnsupportedCommandError,
)


@pytest.mark.unit
def test_build_commands_config_single_command(bcftools_manager, example_sidecar_metadata_inputs_outputs):
    """Test that build_commands_config correctly structures a configuration for a single command."""
    command = "view -s SAMPLES"

    result = bcftools_manager.build_commands_config(command, example_sidecar_metadata_inputs_outputs)

    assert isinstance(result, list)
    assert len(result) == 1, "Should create exactly one command configuration"

    cmd_config = result[0]

    expected_values = {
        "command": command,
        "counter": 0,
        "input_files": example_sidecar_metadata_inputs_outputs["filenames"],
        "sample_subset": example_sidecar_metadata_inputs_outputs["sample_and_filename_subset"],
        "output_temp_files": example_sidecar_metadata_inputs_outputs["output_temp_files"],
    }

    for key, expected_value in expected_values.items():
        if key == "output_temp_files":
            assert len(cmd_config[key]) == len(expected_value)
            for filename in cmd_config[key]:
                assert filename.startswith("temp_subset_") and filename.endswith(".vcf.gz")
        else:
            assert cmd_config[key] == expected_value, (
                f"Expected {key} to be {expected_value}, but got {cmd_config[key]}"
            )

    assert len(cmd_config["output_temp_files"]) == len(example_sidecar_metadata_inputs_outputs["filenames"])


@pytest.mark.unit
def test_build_commands_config_two_commands(bcftools_manager, example_sidecar_metadata_inputs_outputs):
    """Test that build_commands_config correctly structures a configuration for two commands."""
    commands = "view -s SAMPLES; view -r 21:15000000-25000000"

    result = bcftools_manager.build_commands_config(commands, example_sidecar_metadata_inputs_outputs)

    assert isinstance(result, list)
    assert len(result) == 2, "Should create two command configurations"

    expected_configs = [
        {
            "command": "view -s SAMPLES",
            "counter": 0,
            "input_files": example_sidecar_metadata_inputs_outputs["filenames"],
            "sample_subset": example_sidecar_metadata_inputs_outputs["sample_and_filename_subset"],
            "output_temp_files": example_sidecar_metadata_inputs_outputs["output_temp_files"],
            "output_files_count": len(example_sidecar_metadata_inputs_outputs["filenames"]),
        },
        {
            "command": "view -r 21:15000000-25000000",
            "counter": 1,
            "input_files": example_sidecar_metadata_inputs_outputs["output_temp_files"],
            "sample_subset": example_sidecar_metadata_inputs_outputs["sample_and_filename_subset"],
            "output_temp_files": [
                f"temp_subset_1_{i}.vcf.gz"
                for i in range(len(example_sidecar_metadata_inputs_outputs["output_temp_files"]))
            ],
            "output_files_count": len(example_sidecar_metadata_inputs_outputs["output_temp_files"]),
        },
    ]

    for i, cmd_config in enumerate(result):
        expected = expected_configs[i]

        for key in ["command", "counter", "sample_subset"]:
            assert cmd_config[key] == expected[key], (
                f"Command {i + 1}: Expected {key} to be {expected[key]}, but got {cmd_config[key]}"
            )

        for file_key in ["input_files", "output_temp_files"]:
            assert len(cmd_config[file_key]) == len(expected[file_key]), (
                f"Command {i + 1}: Expected {file_key} length to be {len(expected[file_key])}, but got {len(cmd_config[file_key])}"
            )
            for filename, expected_filename in zip(cmd_config[file_key], expected[file_key], strict=False):
                if i == 0 and file_key == "input_files":
                    assert filename == expected_filename, (
                        f"Command {i + 1}: {file_key} file {filename} does not match expected {expected_filename}"
                    )
                else:
                    assert filename.startswith("temp_subset_") and filename.endswith(".vcf.gz"), (
                        f"Command {i + 1}: {file_key} file {filename} does not match expected pattern"
                    )


@pytest.mark.unit
def test_build_commands_config_special_characters(bcftools_manager, example_sidecar_metadata_inputs_outputs):
    """Test that build_commands_config handles edge cases correctly with valid bcftools commands."""

    # Case 1: Command with special characters
    special_char_cmd = "view -s SAMPLES; view --min-af 0.05 --max-af 0.95; view -i 'GT=\"het\"'"
    result = bcftools_manager.build_commands_config(special_char_cmd, example_sidecar_metadata_inputs_outputs)

    assert len(result) == 3, "Should create three command configurations"
    assert result[0]["command"] == "view -s SAMPLES"
    assert result[1]["command"] == "view --min-af 0.05 --max-af 0.95"
    assert result[2]["command"] == "view -i 'GT=\"het\"'"

    # Case 2: Empty command in the middle that should be skipped
    empty_middle_cmd = "view -s SAMPLES;; view -i 'INFO/DP>10'"
    result = bcftools_manager.build_commands_config(empty_middle_cmd, example_sidecar_metadata_inputs_outputs)
    assert len(result) == 2, "Should create two command configurations (empty one is skipped)"
    assert result[0]["command"] == "view -s SAMPLES"
    assert result[1]["command"] == "view -i 'INFO/DP>10'"

    # Case 3: Commands with extra whitespace
    whitespace_cmd = "  view -s SAMPLES  ;  view -G  "
    result = bcftools_manager.build_commands_config(whitespace_cmd, example_sidecar_metadata_inputs_outputs)

    assert len(result) == 2, "Should create two command configurations"
    assert result[0]["command"] == "view -s SAMPLES", "Leading/trailing spaces should be stripped"
    assert result[1]["command"] == "view -G", "Leading/trailing spaces should be stripped"

    # Case 4: Command with quotation marks and complex filtering
    quoted_cmd = "view -s \"SAMPLE1,SAMPLE2\"; view -r 1:1000-2000; view -i 'F_MISSING < 0.1 && MAF[0] > 0.01'"
    result = bcftools_manager.build_commands_config(quoted_cmd, example_sidecar_metadata_inputs_outputs)

    assert len(result) == 3, "Should create three command configurations"
    assert result[0]["command"] == 'view -s "SAMPLE1,SAMPLE2"', "Quotes should be preserved"
    assert result[1]["command"] == "view -r 1:1000-2000"
    assert result[2]["command"] == "view -i 'F_MISSING < 0.1 && MAF[0] > 0.01'"

    # Case 5: Command with multiple consecutive semicolons
    multi_semicolon_cmd = "view -s SAMPLES;;;view -m2 -M2 -v snps"
    result = bcftools_manager.build_commands_config(multi_semicolon_cmd, example_sidecar_metadata_inputs_outputs)
    assert len(result) == 2, "Should create two command configurations (empty ones are skipped)"
    assert result[0]["command"] == "view -s SAMPLES"
    assert result[1]["command"] == "view -m2 -M2 -v snps", "Should select biallelic SNPs"


@pytest.mark.unit
def test_bcftools_pipe_unsupported_command_error(bcftools_manager):
    """Test that BcftoolsPipeUnsupportedCommandError contains the expected attributes and message."""
    command = "merge -m none"
    cmd_name = command.split()[0]
    position = 2
    valid_commands = bcftools_manager.VALID_BCFTOOLS_COMMANDS

    error = BcftoolsPipeUnsupportedCommandError(cmd_name, position, valid_commands)

    assert error.command == cmd_name
    assert error.position == position
    assert error.valid_commands == valid_commands
    expected_message = (
        f"Unsupported bcftools command '{cmd_name}' at position {position}. "
        f"Only the following commands are supported: {', '.join(valid_commands)}"
    )
    assert str(error) == expected_message


@pytest.mark.unit
def test_build_commands_config_invalid_commands(bcftools_manager, example_sidecar_metadata_inputs_outputs):
    """Test that build_commands_config rejects unsupported bcftools commands."""

    merge_command = "view -s SAMPLES; merge -m none; view -i 'GT=\"het\"'"

    with pytest.raises(BcftoolsPipeUnsupportedCommandError) as exc_info:
        bcftools_manager.build_commands_config(merge_command, example_sidecar_metadata_inputs_outputs)

    error_msg = str(exc_info.value)
    assert "Unsupported bcftools command 'merge'" in error_msg
    assert "position 2" in error_msg, "Error should indicate the position of the invalid command"
    assert "merge" in error_msg, "Error should list valid commands"

    multiple_invalid_cmd = "norm --fasta-ref ref.fa; annotate --annotations file.vcf"

    with pytest.raises(BcftoolsPipeUnsupportedCommandError) as exc_info:
        bcftools_manager.build_commands_config(multiple_invalid_cmd, example_sidecar_metadata_inputs_outputs)

    error_msg = str(exc_info.value)
    assert "Unsupported bcftools command 'norm'" in error_msg
    assert "position 1" in error_msg, "Error should indicate the position of the invalid command"


@pytest.mark.unit
def test_build_commands_config_empty_command(bcftools_manager, example_sidecar_metadata_inputs_outputs):
    """Test that fails on various variations of when user inputs empty --command to bcftools-pipe in the CLI."""

    with pytest.raises(BcftoolsPipeEmptyCommandError) as exc_info:
        bcftools_manager.build_commands_config("", example_sidecar_metadata_inputs_outputs)
    assert "Empty command provided" in str(exc_info.value)

    with pytest.raises(BcftoolsPipeEmptyCommandError) as exc_info:
        bcftools_manager.build_commands_config(";", example_sidecar_metadata_inputs_outputs)
    assert "Empty command provided" in str(exc_info.value)

    with pytest.raises(BcftoolsPipeEmptyCommandError) as exc_info:
        bcftools_manager.build_commands_config("  ;  ", example_sidecar_metadata_inputs_outputs)
    assert "Empty command provided" in str(exc_info.value)


@pytest.mark.unit
def test_bcftools_environment_error():
    """Test that BcftoolsEnvironmentError contains the expected attributes and message."""
    container_name = "test-container"
    error = BcftoolsEnvironmentError(container_name)

    assert error.container_name == container_name
    assert f"No running container found with name {container_name}" in str(error)


@pytest.mark.unit
def test_bcftools_command_error():
    """Test that BcftoolsCommandError contains the expected attributes and message."""
    command = "view -h sample.vcf"
    error_details = "Command failed"
    error = BcftoolsCommandError(command=command, error_details=error_details)

    assert error.command == command
    assert error.error_details == error_details

    assert f"bcftools command failed: '{command}'" in str(error)
    assert f"with error details: {error_details}" in str(error)


@pytest.mark.unit
@patch("subprocess.run")
def test_get_container_id_subprocess_error(mock_run, bcftools_manager):
    """Test that get_container_id raises BcftoolsEnvironmentError on subprocess error."""
    mock_run.side_effect = subprocess.SubprocessError("Docker command failed")

    with pytest.raises(BcftoolsEnvironmentError) as excinfo:
        bcftools_manager.get_container_id(bcftools_manager.CONTAINER_NAME)

    assert bcftools_manager.CONTAINER_NAME in str(excinfo.value)

    mock_run.assert_called_with(
        ["docker", "ps", "--filter", f"name={bcftools_manager.CONTAINER_NAME}", "--format", "{{.ID}}"],
        capture_output=True,
        text=True,
        check=True,
    )


@pytest.mark.unit
@patch("divbase_lib.queries.BcftoolsQueryManager.get_container_id")
@patch("os.path.exists")
def test_run_bcftools_container_not_found(mock_exists_in_docker, mock_get_container_id, bcftools_manager):
    """
    Test that BcftoolsEnvironmentError is raised from self.run_bcftools when container is not found.

    The mock_exists_in_docker simulate the synchronous scenario when the code is not executed in a Docker container (i.e., no /.dockerenv file),
    which makes the code look for a valid container ID.

    The `mock_run_get_container_id` simulates the scenario where the command to get the container ID returns an empty string, indicating
    that no container is running with the expected name.

    Together, this should raise a BcftoolsEnvironmentError.
    """
    mock_exists_in_docker.return_value = False
    mock_get_container_id.side_effect = BcftoolsEnvironmentError(bcftools_manager.CONTAINER_NAME)

    with pytest.raises(BcftoolsEnvironmentError) as excinfo:
        bcftools_manager.run_bcftools("view -h sample.vcf")

    assert bcftools_manager.CONTAINER_NAME in str(excinfo.value)

    mock_exists_in_docker.assert_called_with("/.dockerenv")
    mock_get_container_id.assert_called_with(bcftools_manager.CONTAINER_NAME)


@pytest.mark.unit
@patch("subprocess.run")
@patch("divbase_lib.queries.BcftoolsQueryManager.get_container_id")
@patch("os.path.exists")
def test_command_failure_exec_into_container(mock_exists_in_docker, mock_get_container_id, mock_run, bcftools_manager):
    """Test that BcftoolsCommandError is raised when a command fails in the container.

    The mock_exists_in_docker simulate the synchronous scenario when the code is not executed in a Docker container (i.e., no /.dockerenv file),
    which makes the code look for a valid container ID.
    The mock_get_container_id simulates the scenario where the command to get the container ID returns a valid ID, indicating that a container is running.

    Then mock_run.side_effect simulates a failure when running a bcftools command inside the container.
    This raises a CalledProcessError in the test, which is caught and re-raised as a BcftoolsCommandError.
    """

    mock_exists_in_docker.return_value = False
    mock_get_container_id.return_value = "abc123"

    command_error = subprocess.CalledProcessError(returncode=1, cmd=["bcftools", "view", "-h", "sample.vcf"])
    command_error.stderr = "non-zero exit status 1d"
    mock_run.side_effect = command_error

    with pytest.raises(BcftoolsCommandError) as excinfo:
        bcftools_manager.run_bcftools("view -h sample.vcf")

    assert "view -h sample.vcf" in str(excinfo.value)
    assert "non-zero exit status 1" in str(excinfo.value)


@pytest.mark.unit
@patch("subprocess.run")
@patch("os.path.exists")
def test_command_failure_async_inside_container(mock_exists_in_docker, mock_run, bcftools_manager):
    """
    Test that BcftoolsCommandError is raised when a command fails inside a container after being recived from the task queue.

    The mock_exists_in_docker simulate the asynchronous scenario when the task is picked up by a worker and executed inside
    (i.e. /.dockerenv file exists).

    The mock_run then simulates the scenario where the command to run bcftools fails, raising a BcftoolsCommandError.

    """
    mock_exists_in_docker.return_value = True

    command_error = subprocess.CalledProcessError(returncode=1, cmd=["bcftools", "view", "-h", "sample.vcf"])
    command_error.stderr = "non-zero exit status 1"
    mock_run.side_effect = command_error

    with pytest.raises(BcftoolsCommandError) as excinfo:
        bcftools_manager.run_bcftools("view -h sample.vcf")

    assert "view -h sample.vcf" in str(excinfo.value)
    assert "non-zero exit status 1" in str(excinfo.value)


@pytest.mark.unit
def test_temp_file_management_cleans_up_files(bcftools_manager):
    """
    Test that temp_file_management context manager cleans up temporary files.
    To simulate that the `finally` block of the context manager is executed, bcftools_manager.temp_file_management()
    is called with a `with` statement. The cleanup of temporary files (i.e. the `finally` block) should then be
    triggered when the `with` block is exited.
    """
    bcftools_manager.cleanup_temp_files = MagicMock()

    with bcftools_manager.temp_file_management():
        bcftools_manager.temp_files.extend(["file1.tmp", "file2.tmp"])

    bcftools_manager.cleanup_temp_files.assert_called_once_with(["file1.tmp", "file2.tmp"])
