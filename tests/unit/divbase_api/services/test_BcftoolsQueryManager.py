"""Unit tests for BcftoolsQueryManager class, i.e. service layer called by the vcf query task."""

import subprocess
from unittest.mock import MagicMock, patch

import pytest

from divbase_api.services.bcftools_helpers import get_container_id, run_bcftools
from divbase_lib.exceptions import (
    BcftoolsCommandError,
    BcftoolsEnvironmentError,
    BcftoolsPipeEmptyCommandError,
    BcftoolsPipeUnsupportedCommandError,
)


@pytest.mark.unit
def test_build_commands_config_single_command(bcftools_manager, example_sidecar_metadata_inputs_outputs):
    """Test that build_commands_config correctly structures a configuration for a single command."""
    command = "view -r 21:15000000-25000000"

    bcftools_inputs = example_sidecar_metadata_inputs_outputs["bcftools_inputs"]
    result = bcftools_manager.build_commands_config(command, bcftools_inputs)

    assert isinstance(result, list)
    assert len(result) == 1, "Should create exactly one command configuration"

    cmd_config = result[0]

    assert cmd_config.command == command
    assert cmd_config.counter == 0
    assert cmd_config.input_files == bcftools_inputs.filenames
    assert cmd_config.sample_subset == bcftools_inputs.sample_and_filename_subset
    assert len(cmd_config.output_temp_files) == len(example_sidecar_metadata_inputs_outputs["output_temp_files"])
    for filename in cmd_config.output_temp_files:
        assert filename.startswith("temp_subset_") and filename.endswith(".bcf")
    assert len(cmd_config.output_temp_files) == len(bcftools_inputs.filenames)
    assert cmd_config.auto_sample_injection is True


@pytest.mark.unit
def test_build_commands_config_two_commands(bcftools_manager, example_sidecar_metadata_inputs_outputs):
    """Test that build_commands_config correctly structures a configuration for two commands."""
    commands = "view -r 21:15000000-25000000; view -G"

    bcftools_inputs = example_sidecar_metadata_inputs_outputs["bcftools_inputs"]
    result = bcftools_manager.build_commands_config(commands, bcftools_inputs)

    assert isinstance(result, list)
    assert len(result) == 2, "Should create two command configurations"

    cmd1 = result[0]
    assert cmd1.command == "view -r 21:15000000-25000000"
    assert cmd1.counter == 0
    assert cmd1.sample_subset == bcftools_inputs.sample_and_filename_subset
    assert cmd1.input_files == bcftools_inputs.filenames
    assert len(cmd1.output_temp_files) == len(example_sidecar_metadata_inputs_outputs["output_temp_files"])
    for filename in cmd1.output_temp_files:
        assert filename.startswith("temp_subset_") and filename.endswith(".bcf")
    assert cmd1.auto_sample_injection is True

    cmd2 = result[1]
    assert cmd2.command == "view -G"
    assert cmd2.counter == 1
    assert cmd2.sample_subset == bcftools_inputs.sample_and_filename_subset
    assert cmd2.input_files == cmd1.output_temp_files
    assert len(cmd2.output_temp_files) == len(example_sidecar_metadata_inputs_outputs["output_temp_files"])
    for filename in cmd2.output_temp_files:
        assert filename.startswith("temp_subset_") and filename.endswith(".bcf")
    assert cmd2.auto_sample_injection is True


@pytest.mark.unit
def test_build_commands_config_respects_auto_sample_injection_flag(
    bcftools_manager, example_sidecar_metadata_inputs_outputs
):
    """Test that auto sample injection can be disabled (all-samples mode behavior)."""

    original_inputs = example_sidecar_metadata_inputs_outputs["bcftools_inputs"]
    bcftools_inputs = type(original_inputs)(
        filenames=original_inputs.filenames,
        sample_and_filename_subset=original_inputs.sample_and_filename_subset,
        sampleIDs=original_inputs.sampleIDs,
        auto_sample_injection=False,
    )
    result = bcftools_manager.build_commands_config("view -r 21:15000000-25000000", bcftools_inputs)

    assert len(result) == 1
    assert result[0].auto_sample_injection is False


@pytest.mark.unit
def test_build_commands_config_special_characters(bcftools_manager, example_sidecar_metadata_inputs_outputs):
    """Test that build_commands_config handles edge cases correctly with valid bcftools commands."""
    bcftools_inputs = example_sidecar_metadata_inputs_outputs["bcftools_inputs"]

    # Case 1: Command with special characters
    special_char_cmd = "view -r 1:1000-2000; view --min-af 0.05 --max-af 0.95; view -i 'GT=\"het\"'"
    result = bcftools_manager.build_commands_config(special_char_cmd, bcftools_inputs)

    assert len(result) == 3, "Should create three command configurations"
    assert result[0].command == "view -r 1:1000-2000"
    assert result[1].command == "view --min-af 0.05 --max-af 0.95"
    assert result[2].command == "view -i 'GT=\"het\"'"

    # Case 2: Empty command in the middle that should be skipped
    empty_middle_cmd = "view -r 1:1000-2000;; view -i 'INFO/DP>10'"
    result = bcftools_manager.build_commands_config(empty_middle_cmd, bcftools_inputs)
    assert len(result) == 2, "Should create two command configurations (empty one is skipped)"
    assert result[0].command == "view -r 1:1000-2000"
    assert result[1].command == "view -i 'INFO/DP>10'"

    # Case 3: Commands with extra whitespace
    whitespace_cmd = "  view -r 1:1000-2000  ;  view -G  "
    result = bcftools_manager.build_commands_config(whitespace_cmd, bcftools_inputs)

    assert len(result) == 2, "Should create two command configurations"
    assert result[0].command == "view -r 1:1000-2000", "Leading/trailing spaces should be stripped"
    assert result[1].command == "view -G", "Leading/trailing spaces should be stripped"

    # Case 4: Command with quotation marks and complex filtering
    quoted_cmd = "view -r \"1:1000-2000\"; view -r 2:1000-2000; view -i 'F_MISSING < 0.1 && MAF[0] > 0.01'"
    result = bcftools_manager.build_commands_config(quoted_cmd, bcftools_inputs)

    assert len(result) == 3, "Should create three command configurations"
    assert result[0].command == 'view -r "1:1000-2000"', "Quotes should be preserved"
    assert result[1].command == "view -r 2:1000-2000"
    assert result[2].command == "view -i 'F_MISSING < 0.1 && MAF[0] > 0.01'"

    # Case 5: Command with multiple consecutive semicolons
    multi_semicolon_cmd = "view -r 1:1000-2000;;;view -m2 -M2 -v snps"
    result = bcftools_manager.build_commands_config(multi_semicolon_cmd, bcftools_inputs)
    assert len(result) == 2, "Should create two command configurations (empty ones are skipped)"
    assert result[0].command == "view -r 1:1000-2000"
    assert result[1].command == "view -m2 -M2 -v snps", "Should select biallelic SNPs"


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
    bcftools_inputs = example_sidecar_metadata_inputs_outputs["bcftools_inputs"]

    merge_command = "view -r 1:1000-2000; merge -m none; view -i 'GT=\"het\"'"

    with pytest.raises(BcftoolsPipeUnsupportedCommandError) as exc_info:
        bcftools_manager.build_commands_config(merge_command, bcftools_inputs)

    error_msg = str(exc_info.value)
    assert "Unsupported bcftools command 'merge'" in error_msg
    assert "position 2" in error_msg, "Error should indicate the position of the invalid command"
    assert "merge" in error_msg, "Error should list valid commands"

    multiple_invalid_cmd = "norm --fasta-ref ref.fa; annotate --annotations file.vcf"

    with pytest.raises(BcftoolsPipeUnsupportedCommandError) as exc_info:
        bcftools_manager.build_commands_config(multiple_invalid_cmd, bcftools_inputs)

    error_msg = str(exc_info.value)
    assert "Unsupported bcftools command 'norm'" in error_msg
    assert "position 1" in error_msg, "Error should indicate the position of the invalid command"


@pytest.mark.unit
def test_build_commands_config_empty_command(bcftools_manager, example_sidecar_metadata_inputs_outputs):
    """Test that fails on various variations of when user inputs empty --command to divbase-cli query vcf command in the CLI."""
    bcftools_inputs = example_sidecar_metadata_inputs_outputs["bcftools_inputs"]

    with pytest.raises(BcftoolsPipeEmptyCommandError) as exc_info:
        bcftools_manager.build_commands_config("", bcftools_inputs)
    assert "Empty command provided" in str(exc_info.value)

    with pytest.raises(BcftoolsPipeEmptyCommandError) as exc_info:
        bcftools_manager.build_commands_config(";", bcftools_inputs)
    assert "Empty command provided" in str(exc_info.value)

    with pytest.raises(BcftoolsPipeEmptyCommandError) as exc_info:
        bcftools_manager.build_commands_config("  ;  ", bcftools_inputs)
    assert "Empty command provided" in str(exc_info.value)


@pytest.mark.unit
@patch("os.path.exists", return_value=True)
def test_execute_pipe_empty_command_raises_before_processing(
    mock_exists_in_docker, bcftools_manager, example_sidecar_metadata_inputs_outputs
):
    """
    Test that if execute_pipe gets an empty command string, it should fail immediately.
    """
    bcftools_inputs = example_sidecar_metadata_inputs_outputs["bcftools_inputs"]
    with pytest.raises(BcftoolsPipeEmptyCommandError):
        bcftools_manager.execute_pipe("", bcftools_inputs, job_id=1)


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
        get_container_id(bcftools_manager.CONTAINER_NAME)

    assert bcftools_manager.CONTAINER_NAME in str(excinfo.value)

    mock_run.assert_called_with(
        ["docker", "ps", "--filter", f"name={bcftools_manager.CONTAINER_NAME}", "--format", "{{.ID}}"],
        capture_output=True,
        text=True,
        check=True,
    )


@pytest.mark.unit
@patch("divbase_api.services.bcftools_helpers.get_container_id")
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
        run_bcftools("view -h sample.vcf")

    assert bcftools_manager.CONTAINER_NAME in str(excinfo.value)

    mock_exists_in_docker.assert_called_with("/.dockerenv")
    mock_get_container_id.assert_called_with(bcftools_manager.CONTAINER_NAME)


@pytest.mark.unit
@patch("subprocess.Popen")
@patch("divbase_api.services.bcftools_helpers.get_container_id")
@patch("os.path.exists")
def test_command_failure_exec_into_container(
    mock_exists_in_docker, mock_get_container_id, mock_popen, bcftools_manager
):
    """Test that BcftoolsCommandError is raised when a command fails in the container.

    The mock_exists_in_docker simulate the synchronous scenario when the code is not executed in a Docker container (i.e., no /.dockerenv file),
    which makes the code look for a valid container ID.
    The mock_get_container_id simulates the scenario where the command to get the container ID returns a valid ID, indicating that a container is running.

    Then mock_run.side_effect simulates a failure when running a bcftools command inside the container.
    This raises a CalledProcessError in the test, which is caught and re-raised as a BcftoolsCommandError.
    """

    mock_exists_in_docker.return_value = False
    mock_get_container_id.return_value = "abc123"

    mock_popen.side_effect = OSError("non-zero exit status 1")

    with pytest.raises(BcftoolsCommandError) as excinfo:
        run_bcftools("view -h sample.vcf")

    assert "view -h sample.vcf" in str(excinfo.value)
    assert "non-zero exit status 1" in str(excinfo.value)


@pytest.mark.unit
@patch("subprocess.Popen")
@patch("os.path.exists")
def test_command_failure_async_inside_container(mock_exists_in_docker, mock_popen, bcftools_manager):
    """
    Test that BcftoolsCommandError is raised when a command fails inside a container after being recived from the task queue.

    The mock_exists_in_docker simulate the asynchronous scenario when the task is picked up by a worker and executed inside
    (i.e. /.dockerenv file exists).

    The mock_run then simulates the scenario where the command to run bcftools fails, raising a BcftoolsCommandError.

    """
    mock_exists_in_docker.return_value = True

    mock_popen.side_effect = OSError("non-zero exit status 1")

    with pytest.raises(BcftoolsCommandError) as excinfo:
        run_bcftools("view -h sample.vcf")

    assert "view -h sample.vcf" in str(excinfo.value)
    assert "non-zero exit status 1" in str(excinfo.value)


@pytest.mark.unit
@patch("subprocess.Popen")
@patch("os.path.exists")
def test_run_bcftools_uses_shlex_split_for_quoted_args(mock_exists_in_docker, mock_popen, bcftools_manager):
    """Test that quoted expressions passed as single commands/args to bcftools subprocess call."""

    mock_exists_in_docker.return_value = True
    mock_proc = MagicMock()
    mock_popen.return_value = mock_proc

    command = "view -i 'GT=\"het\"' sample.vcf"
    result = run_bcftools(command)

    assert result == mock_proc
    mock_popen.assert_called_once_with(["bcftools", "view", "-i", 'GT="het"', "sample.vcf"])


@pytest.mark.unit
def test_run_bcftools_parse_error_raises_bcftools_command_error(bcftools_manager):
    """Test that malformed shell-like command strings raise BcftoolsCommandError."""

    with pytest.raises(BcftoolsCommandError) as excinfo:
        run_bcftools('view -i \'GT="het" sample.vcf')

    assert "Could not parse bcftools command arguments" in str(excinfo.value)


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
