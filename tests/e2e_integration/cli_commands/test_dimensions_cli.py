"""
Tests for the "divbase-cli dimensions" subcommand
"""

import ast
import gzip
import os
import re
from unittest.mock import MagicMock, patch

import pytest
import yaml
from sqlalchemy import delete, select
from typer.testing import CliRunner

from divbase_api.models.vcf_dimensions import VCFMetadataDB, VCFMetadataSamplesDB, VCFMetadataScaffoldsDB
from divbase_api.services.s3_client import create_s3_file_manager
from divbase_api.worker.crud_dimensions import (
    SkippedVCFData,
    create_or_update_skipped_vcf,
    delete_skipped_vcf_batch,
    delete_vcf_metadata,
    delete_vcf_metadata_batch,
    get_skipped_vcfs_by_project_worker,
    get_vcf_metadata_by_project,
)
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


@patch("divbase_api.worker.tasks.create_s3_file_manager")
def test_update_dimensions_cleans_up_index_when_all_vcfs_deleted_from_bucket(
    mock_create_s3_manager,
    CONSTANTS,
    run_update_dimensions,
    db_session_sync,
    project_map,
):
    """
    Test for the case where:
    1. A user uploads VCFs and runs dimensions update (files get indexed),
    2. Then deletes all VCFs from the bucket,
    3. Then runs dimensions update again.

    This should result in the stale DB entries being deleted and a message indicating the files were deleted being sent in the task results.
    """
    mock_create_s3_manager.side_effect = lambda url=None: create_s3_file_manager(url=CONSTANTS["MINIO_URL"])

    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    project_id = project_map[project_name]
    user_id = 1

    # First run: index VCFs normally using real MinIO
    first_result = run_update_dimensions(
        bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
    )
    assert first_result["status"] == "completed"
    indexed_files = first_result.get("VCF_files_added") or []
    assert indexed_files, "Expected at least one indexed VCF file from the first run"

    # Verify the index has entries in the DB before the second run
    vcf_dimensions_before = get_vcf_metadata_by_project(project_id=project_id, db=db_session_sync)
    assert len(vcf_dimensions_before.get("vcf_files", [])) > 0, "Expected DB to have indexed entries"

    # Second run: simulate all VCFs deleted from the bucket by returning an empty S3
    mock_empty_s3 = MagicMock()
    mock_empty_s3.list_files.return_value = []
    mock_empty_s3.latest_version_of_all_files.return_value = {}
    mock_create_s3_manager.side_effect = lambda url=None: mock_empty_s3

    # Regression test: Assert that NoVCFFilesFoundError is not raised when stale index entries exist
    # (earlier versions of the code raised NoVCFFilesFoundError before cleaning up the index)
    try:
        second_result = update_vcf_dimensions_task(
            bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
        )
    except NoVCFFilesFoundError as e:
        raise AssertionError("NoVCFFilesFoundError should not be raised when cleaning up stale index entries") from e

    assert second_result["status"] == "completed", f"Expected completed status, got: {second_result}"
    assert second_result["VCF_files_added"] is None, "Expected no new files indexed"
    assert second_result["VCF_files_skipped"] is None, "Expected no skipped files"
    deleted_files = second_result.get("VCF_files_deleted") or []
    assert len(deleted_files) > 0, "Expected previously indexed files to be reported as deleted"
    for file in indexed_files:
        assert file in deleted_files, f"Expected {file} to be in VCF_files_deleted: {deleted_files}"

    # DB index should now be empty for this project
    db_session_sync.expire_all()
    vcf_dimensions_after = get_vcf_metadata_by_project(project_id=project_id, db=db_session_sync)
    assert vcf_dimensions_after.get("vcf_files") == [], (
        f"Expected DB index to be empty after cleanup, got: {vcf_dimensions_after.get('vcf_files')}"
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


def test_delete_vcf_metadata_batch(
    CONSTANTS,
    run_update_dimensions,
    db_session_sync,
    project_map,
):
    """
    Test that delete_vcf_metadata_batch removes the specified entries and leaves others intact.
    """
    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    project_id = project_map[project_name]
    user_id = 1

    first_result = run_update_dimensions(
        bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
    )
    indexed_files = first_result.get("VCF_files_added", [])
    assert len(indexed_files) >= 2, "Need at least 2 indexed files for this test"

    batch_to_delete = indexed_files[:2]
    remaining_files = indexed_files[2:]

    delete_vcf_metadata_batch(db=db_session_sync, vcf_file_s3_key_batch=batch_to_delete, project_id=project_id)

    db_session_sync.expire_all()
    updated_dimensions = get_vcf_metadata_by_project(project_id=project_id, db=db_session_sync)
    filenames = [entry["vcf_file_s3_key"] for entry in updated_dimensions.get("vcf_files", [])]

    for vcf_file in batch_to_delete:
        assert vcf_file not in filenames, f"Expected {vcf_file} to be deleted, but found in: {filenames}"

    for vcf_file in remaining_files:
        assert vcf_file in filenames, f"Expected {vcf_file} to remain, but missing from: {filenames}"


def test_delete_skipped_vcf_batch(
    CONSTANTS,
    db_session_sync,
    project_map,
):
    """
    Test that delete_skipped_vcf_batch removes the specified skipped VCF entries and leaves others intact.
    """
    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    project_id = project_map[project_name]

    skipped_files = ["skipped_test_1.vcf.gz", "skipped_test_2.vcf.gz", "skipped_test_3.vcf.gz"]
    for vcf_file in skipped_files:
        create_or_update_skipped_vcf(
            db=db_session_sync,
            skipped_vcf_data=SkippedVCFData(
                vcf_file_s3_key=vcf_file,
                project_id=project_id,
                s3_version_id="test-version-id",
                skip_reason="test skip reason",
            ),
        )

    skipped_before = get_skipped_vcfs_by_project_worker(db=db_session_sync, project_id=project_id)
    for vcf_file in skipped_files:
        assert vcf_file in skipped_before, f"Expected {vcf_file} to be present before batch delete"

    batch_to_delete = skipped_files[:2]
    remaining = skipped_files[2:]

    delete_skipped_vcf_batch(db=db_session_sync, vcf_file_s3_key_batch=batch_to_delete, project_id=project_id)

    skipped_after = get_skipped_vcfs_by_project_worker(db=db_session_sync, project_id=project_id)

    for vcf_file in batch_to_delete:
        assert vcf_file not in skipped_after, f"Expected {vcf_file} to be deleted, but still present"

    for vcf_file in remaining:
        assert vcf_file in skipped_after, f"Expected {vcf_file} to remain, but missing"


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


def test_update_dimensions_reindexes_when_child_rows_missing(
    CONSTANTS,
    run_update_dimensions,
    db_session_sync,
    project_map,
):
    """
    Test that an indexed VCF with missing child rows is detected as incomplete and re-indexed.
    """
    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    project_id = project_map[project_name]
    user_id = 1

    first_result = run_update_dimensions(
        bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
    )
    indexed_files = first_result.get("VCF_files_added", [])
    assert indexed_files, "Expected at least one indexed VCF file"
    target_vcf = indexed_files[0]

    target_stmt = select(VCFMetadataDB).where(
        VCFMetadataDB.project_id == project_id, VCFMetadataDB.vcf_file_s3_key == target_vcf
    )
    target_entry = db_session_sync.execute(target_stmt).scalar_one()

    db_session_sync.execute(delete(VCFMetadataSamplesDB).where(VCFMetadataSamplesDB.vcf_metadata_id == target_entry.id))
    db_session_sync.execute(
        delete(VCFMetadataScaffoldsDB).where(VCFMetadataScaffoldsDB.vcf_metadata_id == target_entry.id)
    )
    db_session_sync.commit()

    dimensions_after_child_row_delete = get_vcf_metadata_by_project(project_id=project_id, db=db_session_sync)
    for entry in dimensions_after_child_row_delete["vcf_files"]:
        if entry["vcf_file_s3_key"] == target_vcf:
            incomplete_entry = entry
            break
    assert incomplete_entry["sample_count"] > 0
    assert incomplete_entry["samples"] == []
    assert incomplete_entry["variant_count"] > 0
    assert incomplete_entry["scaffolds"] == []

    second_result = run_update_dimensions(
        bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
    )
    assert second_result["status"] == "completed"
    assert target_vcf in (second_result.get("VCF_files_added") or []), (
        f"Expected incomplete file to be re-indexed, got: {second_result.get('VCF_files_added')}"
    )

    db_session_sync.expire_all()  # Force refresh of db session to avoid stale cache from the session used by run_update_dimensions()
    dimensions_after_reindex = get_vcf_metadata_by_project(project_id=project_id, db=db_session_sync)
    for entry in dimensions_after_reindex["vcf_files"]:
        if entry["vcf_file_s3_key"] == target_vcf:
            repaired_entry = entry
            break

    assert repaired_entry["sample_count"] > 0
    assert len(repaired_entry["samples"]) > 0
    assert repaired_entry["variant_count"] > 0
    assert len(repaired_entry["scaffolds"]) > 0


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


def test_show_dimensions_sample_names_output_writes_file(
    CONSTANTS,
    run_update_dimensions,
    project_map,
    logged_in_edit_user_with_existing_config,
    tmp_path,
):
    """
    Test that --sample-names-output writes full sample-name rows to file.
    """
    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    project_id = project_map[project_name]
    user_id = 1

    run_update_dimensions(bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id)

    output_path = tmp_path / "sample_names.tsv"
    command = f"dimensions show --project {project_name} --sample-names-output {output_path}"
    cli_result = runner.invoke(app, command)
    assert cli_result.exit_code == 0, f"Command failed with: {cli_result.stdout}"
    assert output_path.exists(), f"Expected output file {output_path} to exist"

    rows = output_path.read_text().strip().splitlines()
    assert len(rows) > 0, "Expected at least one sample row in output file"
    assert all("\t" in row for row in rows), "Expected tab-delimited rows in format: filename<TAB>sample_name"
    assert any("HOM_20ind_17SNPs" in row for row in rows), "Expected known VCF filename in output rows"
    assert "Wrote" in cli_result.stdout and "sample-name rows" in cli_result.stdout


def test_show_dimensions_sample_names_stdout_streams_rows(
    CONSTANTS,
    run_update_dimensions,
    project_map,
    logged_in_edit_user_with_existing_config,
):
    """
    Test that --sample-names-stdout prints filename/sample rows to stdout.
    """
    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    project_id = project_map[project_name]
    user_id = 1

    expected_sample_names = [
        "8_HOM-E57",
        "8_HOM-E59",
        "8_HOM-E64",
        "8_HOM-E74",
        "8_HOM-E78",
        "1a_HOM-G34",
        "5a_HOM-I13",
        "5a_HOM-I14",
        "5a_HOM-I20",
        "5a_HOM-I21",
        "5a_HOM-I7",
        "1b_HOM-G55",
        "1b_HOM-G58",
        "1b_HOM-G83",
        "5b_HOM-H17",
        "5b_HOM-H23",
        "5b_HOM-H25",
        "5b_HOM-H7",
        "7_HOM-J21",
        "4_HOM-P25",
    ]

    run_update_dimensions(bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id)

    command = f"dimensions show --project {project_name} --sample-names-stdout"
    cli_result = runner.invoke(app, command)
    assert cli_result.exit_code == 0, f"Command failed with: {cli_result.stdout}"

    output_sample_names = [
        line.split()[1]
        for line in cli_result.stdout.splitlines()
        if len(line.split()) == 2 and line.split()[0].endswith((".vcf", ".vcf.gz"))
    ]
    assert output_sample_names, "Expected streamed sample rows in stdout"
    missing = set(expected_sample_names) - set(output_sample_names)
    assert not missing, f"Missing sample names in output: {missing}"


def test_show_dimensions_truncates_sample_names_in_terminal(
    CONSTANTS,
    run_update_dimensions,
    project_map,
    logged_in_edit_user_with_existing_config,
):
    """
    Test that --sample-names-limit truncates shown sample names and adds a note.
    """
    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    project_id = project_map[project_name]
    user_id = 1

    run_update_dimensions(bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id)

    command = f"dimensions show --project {project_name} --sample-names-limit 2"
    cli_result = runner.invoke(app, command)
    assert cli_result.exit_code == 0, f"Command failed with: {cli_result.stdout}"

    dimensions_info = yaml.safe_load(cli_result.stdout)
    indexed_files = dimensions_info.get("indexed_files", [])
    assert len(indexed_files) > 0, "Expected indexed files in output"

    first_entry_dimensions = indexed_files[0].get("dimensions", {})
    shown_sample_names = first_entry_dimensions.get("sample_names", [])
    assert len(shown_sample_names) <= 2, f"Expected sample_names to be truncated to <=2, got {shown_sample_names}"
    assert "sample_names_note" in first_entry_dimensions, "Expected truncation note in output"


def test_show_dimensions_rejects_output_and_stdout_together(
    CONSTANTS,
    run_update_dimensions,
    project_map,
    logged_in_edit_user_with_existing_config,
    tmp_path,
):
    """
    Test that using --sample-names-output and --sample-names-stdout together fails.
    """
    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    project_id = project_map[project_name]
    user_id = 1

    run_update_dimensions(bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id)

    output_path = tmp_path / "sample_names.tsv"
    command = f"dimensions show --project {project_name} --sample-names-output {output_path} --sample-names-stdout"
    cli_result = runner.invoke(app, command)
    assert cli_result.exit_code != 0, "Expected command to fail when both output modes are provided"
    combined_output = (
        (cli_result.stdout or "")
        + (getattr(cli_result, "stderr", "") or "")
        + (str(cli_result.exception) if cli_result.exception else "")
    )
    assert "Use only one of --sample-names-output or --sample-names-stdout." in combined_output


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
    assert re.search(r"Total columns:\s+\d+", cli_result.stdout), "Expected total columns in summary"
    assert re.search(r"Samples matching project VCF dimensions:\s+\d+/\d+", cli_result.stdout), (
        "Expected dimensions sample match counts in summary"
    )
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
    Only assert fixture-driven messages. Ignore dimensions mismatch/project state-dependent messages.
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
    assert "ERRORS" in cli_result.stdout, "Expected errors section"
    assert "WARNINGS" in cli_result.stdout, "Expected warnings section"

    assert "Expected 4 tab-separated columns from reading the header, found 2" in cli_result.stdout
    # Rich/terminal wrapping can split long phrases across lines, so assert key fragments.
    assert "mixed element types" in cli_result.stdout and "in lists" in cli_result.stdout
    assert "Found 2 cell(s) with invalid list syntax or not parsed as list" in cli_result.stdout
    assert "Sample_ID is empty or missing in 2 row(s)" in cli_result.stdout
    assert "Duplicate Sample_IDs found: 'test_duplicate' appears in 2 row(s)" in cli_result.stdout
    assert "Sample_ID column contains list values" in cli_result.stdout
    assert "Found 3 cell(s) with leading or trailing whitespace" in cli_result.stdout
    assert "This column contains mixed-type cells" in cli_result.stdout
    assert "Validation failed! Please fix the errors above before uploading." in cli_result.stdout


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

    assert cli_result.exit_code == 2, "Expected exit code 2 for nonexistent file (Typer path validation)"
    assert "does not exist" in cli_result.output.lower(), "Expected error message about file not existing"


@patch("divbase_api.worker.tasks.create_s3_file_manager")
def test_update_dimensions_cleans_up_csi_index_files_from_worker(
    mock_create_s3_manager,
    CONSTANTS,
    project_map,
    tmp_path,
    monkeypatch,
):
    """
    Test that CSI index files created during dimension calculation are deleted from the worker
    filesystem after update_vcf_dimensions_task completes.
    """
    mock_create_s3_manager.side_effect = lambda url=None: create_s3_file_manager(url=CONSTANTS["MINIO_URL"])
    monkeypatch.chdir(tmp_path)

    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    project_id = project_map[project_name]
    user_id = 1

    result = update_vcf_dimensions_task(
        bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
    )

    assert result["status"] == "completed"

    vcf_files_remaining = list(tmp_path.glob("*.vcf")) + list(tmp_path.glob("*.vcf.gz"))
    assert vcf_files_remaining == [], (
        f"Expected all VCF files to be deleted from worker after task, but found: {vcf_files_remaining}"
    )

    csi_files_remaining = list(tmp_path.glob("*.csi"))
    assert csi_files_remaining == [], (
        f"Expected all CSI index files to be deleted from worker after task, but found: {csi_files_remaining}"
    )


@patch("divbase_api.worker.tasks.create_s3_file_manager")
def test_update_dimensions_indexes_uncompressed_vcf(
    mock_create_s3_manager,
    CONSTANTS,
    db_session_sync,
    project_map,
    tmp_path,
):
    """
    Test that update_vcf_dimensions_task correctly indexes a plain uncompressed .vcf file.

    bcftools index --csi requires bgzipped input, so calculate_dimensions bgzips plain .vcf
    files to a temp .vcf.gz internally before indexing, then cleans up the temp files.
    This tests that the full indexing pipeline works end-to-end for uncompressed VCFs.
    """
    mock_create_s3_manager.side_effect = lambda url=None: create_s3_file_manager(url=CONSTANTS["MINIO_URL"])

    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    project_id = project_map[project_name]
    user_id = 1

    plain_vcf_name = "test_uncompressed_plain.vcf"
    vcf_path = tmp_path / plain_vcf_name
    vcf_content = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=1,length=22053058>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\n"
        "1\t17504018\t.\tC\tA\t.\tPASS\t.\tGT\t0/1\t0/0\n"
        "1\t22053057\t.\tG\tA\t.\tPASS\t.\tGT\t1/1\t0/1\n"
    )
    vcf_path.write_text(vcf_content)

    s3_file_manager = create_s3_file_manager(url=CONSTANTS["MINIO_URL"])
    s3_file_manager.upload_files(
        to_upload={plain_vcf_name: vcf_path},
        bucket_name=bucket_name,
    )

    try:
        result = update_vcf_dimensions_task(
            bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id
        )

        assert result["status"] == "completed"
        indexed_files = result.get("VCF_files_added") or []
        assert plain_vcf_name in indexed_files, f"Expected plain .vcf file to be indexed, got: {indexed_files}"

        # Verify the dimensions were stored correctly in the DB
        db_session_sync.expire_all()
        vcf_dimensions = get_vcf_metadata_by_project(project_id=project_id, db=db_session_sync)
        all_indexed_keys = [entry["vcf_file_s3_key"] for entry in vcf_dimensions.get("vcf_files", [])]
        assert plain_vcf_name in all_indexed_keys

        entry = next(e for e in vcf_dimensions["vcf_files"] if e["vcf_file_s3_key"] == plain_vcf_name)
        assert entry["sample_count"] == 2
        assert entry["variant_count"] == 2
        assert "1" in entry["scaffolds"]
        assert "sample1" in entry["samples"]
        assert "sample2" in entry["samples"]
    finally:
        # Remove the uploaded file from the bucket to avoid polluting subsequent tests.
        # auto_clean_dimensions_entries_for_all_projects only cleans the DB, not the bucket.
        s3_file_manager.soft_delete_objects(objects=[plain_vcf_name], bucket_name=bucket_name)
