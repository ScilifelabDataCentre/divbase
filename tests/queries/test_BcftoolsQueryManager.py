import pytest


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
        assert cmd_config[key] == expected_value, f"Expected {key} to be {expected_value}, but got {cmd_config[key]}"

    assert len(cmd_config["output_temp_files"]) == len(example_sidecar_metadata_inputs_outputs["filenames"])


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

        for key in ["command", "counter", "input_files", "sample_subset"]:
            assert cmd_config[key] == expected[key], (
                f"Command {i + 1}: Expected {key} to be {expected[key]}, but got {cmd_config[key]}"
            )

        assert cmd_config["output_temp_files"] == expected["output_temp_files"], (
            f"Command {i + 1}: Expected output_temp_files to be {expected['output_temp_files']}, "
            f"but got {cmd_config['output_temp_files']}"
        )

        assert len(cmd_config["output_temp_files"]) == expected["output_files_count"]


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

    # Case 4: Single semicolon (i.e. two empty commands)
    semicolon_cmd = ";"
    result = bcftools_manager.build_commands_config(semicolon_cmd, example_sidecar_metadata_inputs_outputs)
    assert len(result) == 0, "Should create no command configurations (all empty)"

    # Case 5: Command with quotation marks and complex filtering
    quoted_cmd = "view -s \"SAMPLE1,SAMPLE2\"; view -r 1:1000-2000; view -i 'F_MISSING < 0.1 && MAF[0] > 0.01'"
    result = bcftools_manager.build_commands_config(quoted_cmd, example_sidecar_metadata_inputs_outputs)

    assert len(result) == 3, "Should create three command configurations"
    assert result[0]["command"] == 'view -s "SAMPLE1,SAMPLE2"', "Quotes should be preserved"
    assert result[1]["command"] == "view -r 1:1000-2000"
    assert result[2]["command"] == "view -i 'F_MISSING < 0.1 && MAF[0] > 0.01'"

    # Case 6: Command with multiple consecutive semicolons
    multi_semicolon_cmd = "view -s SAMPLES;;;view -m2 -M2 -v snps"
    result = bcftools_manager.build_commands_config(multi_semicolon_cmd, example_sidecar_metadata_inputs_outputs)
    assert len(result) == 2, "Should create two command configurations (empty ones are skipped)"
    assert result[0]["command"] == "view -s SAMPLES"
    assert result[1]["command"] == "view -m2 -M2 -v snps", "Should select biallelic SNPs"


def test_build_commands_config_invalid_commands(bcftools_manager, example_sidecar_metadata_inputs_outputs):
    """Test that build_commands_config rejects unsupported bcftools commands."""

    merge_command = "view -s SAMPLES; merge -m none; view -i 'GT=\"het\"'"

    with pytest.raises(ValueError) as exc_info:
        bcftools_manager.build_commands_config(merge_command, example_sidecar_metadata_inputs_outputs)

    error_msg = str(exc_info.value)
    assert "Unsupported bcftools command 'merge'" in error_msg
    assert "position 2" in error_msg, "Error should indicate the position of the invalid command"
    assert "view" in error_msg, "Error should list valid commands"

    multiple_invalid_cmd = "norm --fasta-ref ref.fa; annotate --annotations file.vcf"

    with pytest.raises(ValueError) as exc_info:
        bcftools_manager.build_commands_config(multiple_invalid_cmd, example_sidecar_metadata_inputs_outputs)

    error_msg = str(exc_info.value)
    assert "Unsupported bcftools command 'norm'" in error_msg
    assert "position 1" in error_msg, "Error should indicate the position of the invalid command"
