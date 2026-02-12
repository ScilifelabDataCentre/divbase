"""
Tests for the "divbase-cli dimensions" subcommand
"""

import ast
import gzip
import os
import re
from unittest.mock import patch

import pytest
import yaml
from typer.testing import CliRunner

from divbase_api.services.s3_client import create_s3_file_manager
from divbase_api.worker.crud_dimensions import delete_vcf_metadata, get_vcf_metadata_by_project
from divbase_api.worker.tasks import update_vcf_dimensions_task
from divbase_cli.cli_exceptions import DivBaseAPIError
from divbase_cli.divbase_cli import app
from divbase_lib.exceptions import NoVCFFilesFoundError

runner = CliRunner()


@pytest.fixture(autouse=True, scope="function")
def auto_clean_dimensions_entries_for_all_projects(clean_all_projects_dimensions):
    """Enable auto-cleanup of dimensions entries for all tests in this test file."""
    yield


def _parse_list_from_cli_output(stdout: str) -> list:
    """
    Helper function to parse a Python list from CLI output that may span multiple lines.
    """
    lines = stdout.splitlines()
    list_text = ""
    collecting = False

    for line in lines:
        if "[" in line and "count:" not in line:
            collecting = True
        if collecting:
            list_text += line
            if "]" in line:
                break

    assert list_text, f"List not found in output:\n{stdout}"

    list_start = list_text.find("[")
    list_end = list_text.rfind("]") + 1
    return ast.literal_eval(list_text[list_start:list_end])


def test_update_vcf_dimensions_task_directly(
    CONSTANTS,
    run_update_dimensions,
    db_session_sync,
    project_map,
):
    """
    Test that runs the update task and verifies all VCF files are indexed via the API.
    """
    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    project_id = project_map[project_name]
    user_id = 1
    result = run_update_dimensions(
        bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
    )

    vcf_files = [f for f in CONSTANTS["PROJECT_CONTENTS"][project_name] if f.endswith(".vcf.gz") or f.endswith(".vcf")]
    indexed_files = result.get("VCF_files_added", [])

    for vcf_file in vcf_files:
        assert vcf_file in indexed_files, f"{vcf_file} not found in indexed files: {indexed_files}"


def test_show_vcf_dimensions_task(
    CONSTANTS,
    run_update_dimensions,
    db_session_sync,
    project_map,
    logged_in_edit_user_with_existing_config,
):
    """
    Test the CLI show command after indexing dimensions via the API.
    """
    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    project_id = project_map[project_name]
    user_id = 1

    run_update_dimensions(bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id)

    # Basic version of command
    command = f"dimensions show --project {project_name}"
    cli_result = runner.invoke(app, command)
    assert cli_result.exit_code == 0

    dimensions_info = yaml.safe_load(cli_result.stdout)
    assert isinstance(dimensions_info, dict), f"Expected dict, got: {type(dimensions_info)}"
    indexed_files = dimensions_info.get("indexed_files", [])

    vcf_files = [f for f in CONSTANTS["PROJECT_CONTENTS"][project_name] if f.endswith(".vcf.gz") or f.endswith(".vcf")]
    found_files = [entry.get("filename") for entry in indexed_files]
    for vcf_file in vcf_files:
        assert vcf_file in found_files, f"{vcf_file} not found in CLI output:\n{cli_result.stdout}"

    # Unique-scaffolds version of command
    command = f"dimensions show --project {project_name} --unique-scaffolds"
    cli_result = runner.invoke(app, command)
    assert cli_result.exit_code == 0

    lines = cli_result.stdout.splitlines()
    scaffold_line = next((line for line in lines if line.startswith("['")), None)
    assert scaffold_line is not None, f"Scaffold list not found in output:\n{cli_result.stdout}"

    scaffold_names = ast.literal_eval(scaffold_line)
    expected_scaffolds = ["1", "4", "5", "6", "7", "8", "13", "18", "20", "21", "22", "24"]
    assert scaffold_names == expected_scaffolds, f"Expected {expected_scaffolds}, got {scaffold_names}"

    # Filename version of command
    for vcf_file in vcf_files:
        command = f"dimensions show --project {project_name} --filename {vcf_file}"
        cli_result = runner.invoke(app, command)
        assert cli_result.exit_code == 0
        entry = yaml.safe_load(cli_result.stdout)
        assert entry.get("filename") == vcf_file
        match = re.search(r"HOM_20ind_17SNPs\.(\d+)\.vcf\.gz", vcf_file)
        if match:
            scaffold_name = match.group(1)
            scaffolds = entry.get("dimensions", {}).get("scaffolds", [])
            assert scaffold_name in scaffolds, f"{scaffold_name} not found in scaffolds for {vcf_file}: {scaffolds}"


def test_show_vcf_dimensions_task_when_file_missing(
    CONSTANTS,
    db_session_sync,
    project_map,
    logged_in_edit_user_with_existing_config,
):
    """
    Test that the CLI handles the case when no dimensions are indexed (empty database).
    """
    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]

    command = f"dimensions show --project {project_name}"

    result = runner.invoke(app, command)
    assert result.exit_code != 0
    assert isinstance(result.exception, DivBaseAPIError)
    assert project_name in str(result.exception)
    assert "vcf_dimensions_entry_missing_error" in str(result.exception)


def test_get_dimensions_info_returns_empty(
    CONSTANTS,
    db_session_sync,
    project_map,
):
    """
    Test that get_dimensions_info returns empty when no dimensions are indexed in the database.
    """
    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    project_id = project_map[project_name]

    result = get_vcf_metadata_by_project(project_id=project_id, db=db_session_sync)
    assert result["vcf_files"] == []


def test_update_vcf_dimensions_task_raises_no_vcf_files_error(
    CONSTANTS,
    db_session_sync,
    project_map,
):
    """
    Test that the update task raises an error when the bucket has no VCF files.
    """
    test_minio_url = CONSTANTS["MINIO_URL"]
    project_name = "empty-project"
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    project_id = project_map[project_name]
    user_id = 1
    with patch("divbase_api.worker.tasks.create_s3_file_manager") as mock_create_s3_manager:
        mock_create_s3_manager.side_effect = lambda url=None: create_s3_file_manager(url=test_minio_url)
        with pytest.raises(NoVCFFilesFoundError):
            update_vcf_dimensions_task(
                bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
            )


def test_remove_VCF_and_update_dimension_entry(
    CONSTANTS,
    db_session_sync,
    project_map,
):
    """
    Test removing a VCF metadata entry
    """
    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    project_id = project_map[project_name]
    vcf_file = "HOM_20ind_17SNPs.8.vcf.gz"

    delete_vcf_metadata(db=db_session_sync, vcf_file_s3_key=vcf_file, project_id=project_id)
    updated_dimensions = get_vcf_metadata_by_project(project_id=project_id, db=db_session_sync)
    filenames = [entry["vcf_file_s3_key"] for entry in updated_dimensions.get("vcf_files", [])]
    assert vcf_file not in filenames


@patch("divbase_api.worker.tasks.create_s3_file_manager")
def test_update_dimensions_skips_divbase_generated_vcf(
    mock_create_s3_manager,
    CONSTANTS,
    db_session_sync,
    project_map,
    tmp_path,
):
    """
    Test that after running a query (which generates a DivBase result VCF),
    update_vcf_dimensions_task skips that file and returns a skip message.
    """
    mock_create_s3_manager.side_effect = lambda url=None: create_s3_file_manager(url=CONSTANTS["MINIO_URL"])

    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    project_id = project_map[project_name]
    user_id = 1

    divbase_vcf_name = "merged_test_divbase_result.vcf.gz"
    vcf_path = tmp_path / divbase_vcf_name
    divbase_vcf_content = (
        "##fileformat=VCFv4.2\n"
        '##DivBase_created="This is a results file created by a DivBase query; Date=Mon Oct 20 12:00:00 2025"\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n"
        "1\t1000\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n"
    )

    with gzip.open(vcf_path, "wt") as f:
        f.write(divbase_vcf_content)

    s3_file_manager = create_s3_file_manager(url=CONSTANTS["MINIO_URL"])

    s3_file_manager.upload_files(
        to_upload={divbase_vcf_name: vcf_path},
        bucket_name=bucket_name,
    )

    result = update_vcf_dimensions_task(
        bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
    )

    skipped_files = result.get("VCF_files_skipped", [])
    assert any(divbase_vcf_name in msg for msg in skipped_files), (
        f"Expected that this file was skipped; {divbase_vcf_name}; got: {skipped_files}"
    )

    assert result["status"] == "completed"


@patch("divbase_api.worker.tasks.create_s3_file_manager")
def test_update_dimensions_twice_with_no_new_VCF_added_inbetween(
    mock_create_s3_manager,
    CONSTANTS,
    db_session_sync,
    project_map,
    tmp_path,
):
    """
    Test that after running update_vcf_dimensions_task twice with no new VCF files added in between,
    the task returns a message indicating no new files were found.
    """
    mock_create_s3_manager.side_effect = lambda url=None: create_s3_file_manager(url=CONSTANTS["MINIO_URL"])

    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    project_id = project_map[project_name]
    user_id = 1

    result_first_run = update_vcf_dimensions_task(
        bucket_name=bucket_name,
        project_id=project_id,
        project_name=project_name,
        user_id=user_id,
    )

    assert result_first_run["status"] == "completed"
    added_files = result_first_run.get("VCF_files_added", [])
    expected_files = CONSTANTS["PROJECT_CONTENTS"]["split-scaffold-project"]
    expected_vcfs = [f for f in expected_files if f.endswith(".vcf.gz")]
    for vcf in expected_vcfs:
        assert vcf in added_files, f"{vcf} not found in indexed files: {added_files}"

    result_second_run = update_vcf_dimensions_task(
        bucket_name=bucket_name,
        project_id=project_id,
        project_name=project_name,
        user_id=user_id,
    )
    assert result_second_run["status"] == "completed"
    assert result_second_run.get("VCF_files_added") is None or result_second_run.get("VCF_files_added") == [], (
        f"Expected no new files indexed, got: {result_second_run.get('VCF_files_added')}"
    )


def test_show_unique_samples(
    CONSTANTS,
    run_update_dimensions,
    db_session_sync,
    project_map,
    logged_in_edit_user_with_existing_config,
):
    """
    Test the CLI 'dimensions show --unique-samples' command.
    """
    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    project_id = project_map[project_name]
    user_id = 1

    run_update_dimensions(bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id)

    command = f"dimensions show --project {project_name} --unique-samples"
    cli_result = runner.invoke(app, command)
    assert cli_result.exit_code == 0, f"Command failed with: {cli_result.stdout}"

    assert "count:" in cli_result.stdout, "Expected count to be displayed in output"
    assert "Unique sample names found" in cli_result.stdout, "Expected header message"
    assert "[" in cli_result.stdout and "]" in cli_result.stdout, "Expected list output"

    sample_names = _parse_list_from_cli_output(cli_result.stdout)

    assert isinstance(sample_names, list), f"Expected list, got {type(sample_names)}"
    assert len(sample_names) > 0, "Expected at least one sample"

    assert sample_names == sorted(sample_names), f"Samples should be sorted: {sample_names}"


def test_show_unique_scaffolds_dedicated_endpoint(
    CONSTANTS,
    run_update_dimensions,
    db_session_sync,
    project_map,
    logged_in_edit_user_with_existing_config,
):
    """
    Test the CLI 'dimensions show --unique-scaffolds' command using the dedicated endpoint.
    This tests both the CRUD function and the CLI integration.
    """
    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    project_id = project_map[project_name]
    user_id = 1

    run_update_dimensions(bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id)

    command = f"dimensions show --project {project_name} --unique-scaffolds"
    cli_result = runner.invoke(app, command)
    assert cli_result.exit_code == 0, f"Command failed with: {cli_result.stdout}"
    assert "count:" in cli_result.stdout, "Expected count to be displayed in output"

    scaffold_names = _parse_list_from_cli_output(cli_result.stdout)

    expected_scaffolds = ["1", "4", "5", "6", "7", "8", "13", "18", "20", "21", "22", "24"]
    assert scaffold_names == expected_scaffolds, f"Expected {expected_scaffolds}, got {scaffold_names}"

    # Verify numeric scaffolds come first, sorted numerically
    numeric_scaffolds = [s for s in scaffold_names if s.isdigit()]
    assert numeric_scaffolds == sorted(numeric_scaffolds, key=int), "Numeric scaffolds should be sorted numerically"


@pytest.mark.parametrize(
    "option_flag,expected_message,expected_items,verify_sorting",
    [
        (
            "--unique-samples",
            "Unique sample names found",
            [
                "1a_HOM-G34",
                "1b_HOM-G55",
                "1b_HOM-G58",
                "1b_HOM-G83",
                "4_HOM-P25",
                "5a_HOM-I13",
                "5a_HOM-I14",
                "5a_HOM-I20",
                "5a_HOM-I21",
                "5a_HOM-I7",
                "5b_HOM-H17",
                "5b_HOM-H23",
                "5b_HOM-H25",
                "5b_HOM-H7",
                "7_HOM-J21",
                "8_HOM-E57",
                "8_HOM-E59",
                "8_HOM-E64",
                "8_HOM-E74",
                "8_HOM-E78",
            ],
            True,  # Should be sorted alphabetically
        ),
        (
            "--unique-scaffolds",
            "Unique scaffold names found",
            ["1", "4", "5", "6", "7", "8", "13", "18", "20", "21", "22", "24"],
            True,  # Should be sorted numerically then alphabetically
        ),
    ],
)
def test_show_unique_items_parametrized(
    CONSTANTS,
    run_update_dimensions,
    db_session_sync,
    project_map,
    logged_in_edit_user_with_existing_config,
    option_flag,
    expected_message,
    expected_items,
    verify_sorting,
):
    """
    Parametrized test for --unique-samples and --unique-scaffolds options.
    Tests both the CRUD functions and CLI integration.
    """
    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    project_id = project_map[project_name]
    user_id = 1

    run_update_dimensions(bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id)

    command = f"dimensions show --project {project_name} {option_flag}"
    cli_result = runner.invoke(app, command)

    assert cli_result.exit_code == 0, f"Command failed with: {cli_result.stdout}"
    assert "count:" in cli_result.stdout, "Expected count to be displayed in output"
    assert expected_message in cli_result.stdout, f"Expected message '{expected_message}' in output"
    assert "[" in cli_result.stdout and "]" in cli_result.stdout, "Expected list output"

    items = _parse_list_from_cli_output(cli_result.stdout)

    assert isinstance(items, list), f"Expected list, got {type(items)}"
    assert len(items) > 0, f"Expected at least one item in {option_flag} output"

    if expected_items is not None:
        assert items == expected_items, f"Expected {expected_items}, got {items}"

    if verify_sorting and option_flag == "--unique-samples":
        assert items == sorted(items), f"Samples should be sorted alphabetically: {items}"
    elif verify_sorting and option_flag == "--unique-scaffolds":
        numeric_items = [s for s in items if s.isdigit()]
        assert numeric_items == sorted(numeric_items, key=int), "Numeric scaffolds should be sorted numerically"


def test_create_metadata_template(
    CONSTANTS,
    run_update_dimensions,
    db_session_sync,
    project_map,
    logged_in_edit_user_with_existing_config,
    tmp_path,
):
    """
    Test the CLI 'dimensions create-metadata-template' command.
    """
    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    project_id = project_map[project_name]
    user_id = 1

    run_update_dimensions(bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id)
    output_filename = f"test_metadata_{project_name}.tsv"
    output_path = tmp_path / output_filename

    command = f"dimensions create-metadata-template --project {project_name} --output {output_path}"
    cli_result = runner.invoke(app, command)
    assert cli_result.exit_code == 0, f"Command failed with: {cli_result.stdout}"
    assert output_path.exists(), f"Expected output file {output_path} to exist"

    with open(output_path, "r") as f:
        lines = f.readlines()
    assert lines[0].strip() == "#Sample_ID", f"Expected header '#Sample_ID', got {lines[0].strip()}"
    assert len(lines) > 1, "Expected at least one sample in the template"
    stdout_lower = cli_result.stdout.lower()
    assert "unique samples" in stdout_lower or "samples found" in stdout_lower, (
        f"Expected message about unique samples, got: {cli_result.stdout}"
    )
    assert str(output_path) in cli_result.stdout or "written" in stdout_lower, (
        "Expected output filename or confirmation message"
    )


def test_create_metadata_template_with_overwrite_prompt(
    CONSTANTS,
    run_update_dimensions,
    db_session_sync,
    project_map,
    logged_in_edit_user_with_existing_config,
    tmp_path,
):
    """
    Test that create-metadata-template prompts when file exists.
    """
    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    project_id = project_map[project_name]
    user_id = 1

    run_update_dimensions(bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id)

    output_filename = f"test_metadata_{project_name}.tsv"
    output_path = tmp_path / output_filename

    # Create template file (does not exist since before)
    command = f"dimensions create-metadata-template --project {project_name} --output {output_path}"
    cli_result = runner.invoke(app, command)
    assert cli_result.exit_code == 0, f"First creation failed: {cli_result.stdout}"

    # Try to create template file again and decline overwrite
    cli_result = runner.invoke(app, command, input="n\n")
    assert "already exists" in cli_result.stdout, "Expected overwrite prompt"
    assert "not written" in cli_result.stdout.lower() or cli_result.exit_code != 0, (
        "Expected message about file not written or non-zero exit"
    )

    # Try to create template file again and accept overwrite
    cli_result = runner.invoke(app, command, input="y\n")
    assert cli_result.exit_code == 0, f"Expected exit code 0, got {cli_result.exit_code}. Output: {cli_result.stdout}"
    assert "already exists" in cli_result.stdout, "Expected overwrite prompt"


def test_validate_metadata_file_valid(
    CONSTANTS,
    run_update_dimensions,
    db_session_sync,
    project_map,
    logged_in_edit_user_with_existing_config,
):
    """
    Test the CLI 'dimensions validate-metadata-file' command with a valid TSV file.
    """
    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    project_id = project_map[project_name]
    user_id = 1

    run_update_dimensions(bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id)

    fixture_path = os.path.join(
        os.path.dirname(__file__), "../..", "fixtures", "sample_metadata_HOM_chr_split_version.tsv"
    )

    command = f"dimensions validate-metadata-file {fixture_path} --project {project_name}"
    cli_result = runner.invoke(app, command)

    assert cli_result.exit_code == 0, f"Expected validation to succeed with exit code 0, got {cli_result.exit_code}"
    assert "VALIDATION SUMMARY" in cli_result.stdout, "Expected validation summary"
    assert "ERRORS" not in cli_result.stdout, f"Did not expect errors, got: {cli_result.stdout}"


def test_validate_metadata_file_with_errors(
    CONSTANTS,
    run_update_dimensions,
    db_session_sync,
    project_map,
    logged_in_edit_user_with_existing_config,
):
    """
    Test the CLI 'dimensions validate-metadata-file' command with an invalid TSV file.
    """
    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    project_id = project_map[project_name]
    user_id = 1

    run_update_dimensions(bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id)

    fixture_path = os.path.join(
        os.path.dirname(__file__), "../..", "fixtures", "sample_metadata_incorrect_formatting_to_test_tsv_validator.tsv"
    )

    command = f"dimensions validate-metadata-file {fixture_path} --project {project_name}"
    cli_result = runner.invoke(app, command)

    assert cli_result.exit_code == 1, f"Expected validation to fail but it passed: {cli_result.stdout}"
    assert "VALIDATION SUMMARY" in cli_result.stdout, "Expected validation summary"
    assert "ERRORS" in cli_result.stdout or "WARNINGS" in cli_result.stdout, "Expected errors or warnings"


def test_validate_metadata_file_nonexistent(
    CONSTANTS,
    logged_in_edit_user_with_existing_config,
):
    """
    Test that validate-metadata-file handles nonexistent files gracefully.
    """
    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]

    command = f"dimensions validate-metadata-file nonexistent_file.tsv --project {project_name}"
    cli_result = runner.invoke(app, command)

    assert cli_result.exit_code == 1, "Expected exit code 1 for nonexistent file"
    assert "not found" in cli_result.stdout.lower(), "Expected error message about file not found"
