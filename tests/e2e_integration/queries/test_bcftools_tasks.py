import pytest
from celery import current_app
from kombu.connection import Connection

from divbase_api.worker.tasks import (
    _calculate_pairwise_overlap_types_for_sample_sets,
    _check_if_samples_can_be_combined_with_bcftools,
    bcftools_pipe_task,
)


@pytest.fixture(autouse=True, scope="function")
def auto_clean_dimensions_entries_for_all_projects(clean_all_projects_dimensions):
    """Enable auto-cleanup of dimensions entries for all tests in this test file."""
    yield


@pytest.mark.integration
def test_bcftools_pipe_task_with_real_worker(
    wait_for_celery_task_completion,
    bcftools_pipe_kwargs_fixture,
    run_update_dimensions,
    db_session_sync,
    project_map,
    CONSTANTS,
):
    """
    Integration test in which bcftools_pipe_task is run with a real Celery worker.
    Runs locally using the docker compose testing stack defined and performs a real sidecar and bcftools
    query by loading VCF files from the tests/fixtures dir. It was designed for having RabbitMQ as the broker, PostgreSQL as the backend,
    and a custom Celery worker image that has bcftools installed.
    (this test does not download fixture from bucket, since that is handled by the CLI layer)

    Does not assert file download and uploads, since that is handled by tests in tests/cli_commands/test_query_cli.py.
    (but the current version of the task does I/O with files from the bucket, so it is tested indirectly)

    """

    project_name = CONSTANTS["QUERY_PROJECT"]
    project_id = project_map[project_name]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    run_update_dimensions(bucket_name=bucket_name, project_id=project_id, project_name=project_name)
    bcftools_pipe_kwargs_fixture["project_id"] = project_id

    broker_url = current_app.conf.broker_url
    with Connection(broker_url) as conn:
        conn.ensure_connection(max_retries=1)

    async_result = bcftools_pipe_task.apply_async(kwargs=bcftools_pipe_kwargs_fixture)
    task_id = async_result.id
    task_result = wait_for_celery_task_completion(task_id=task_id, max_wait=30)

    user_name = bcftools_pipe_kwargs_fixture.get("user_name")
    assert task_result["status"] == "completed"
    assert task_result["submitter"] == user_name
    assert task_result["output_file"].startswith("merged_") and task_result["output_file"].endswith(".vcf.gz")


@pytest.mark.parametrize(
    "sample_sets_dict,expected_identical_diff_order,expected_partly,expected_non_overlap",
    [
        # Case 1: partly overlapping and non-overlapping
        (
            {
                tuple(["S1", "S2", "S3"]): [],
                tuple(["S3", "S4"]): [],
                tuple(["S5", "S6"]): [],
            },
            [],
            [(tuple(["S1", "S2", "S3"]), tuple(["S3", "S4"]))],
            [(tuple(["S1", "S2", "S3"]), tuple(["S5", "S6"]))],
        ),
        # Case 2: all non-overlapping
        (
            {
                tuple(["A"]): [],
                tuple(["B"]): [],
                tuple(["C"]): [],
            },
            [],
            [],
            [
                (tuple(["A"]), tuple(["B"])),
                (tuple(["A"]), tuple(["C"])),
                (tuple(["B"]), tuple(["C"])),
            ],
        ),
        # Case 3: identical elements, different order
        (
            {
                tuple(["X", "Y"]): [],
                tuple(["Y", "X"]): [],
            },
            [(tuple(["X", "Y"]), tuple(["Y", "X"]))],
            [],
            [],
        ),
    ],
)
def test_calculate_pairwise_overlap_types_for_sample_sets(
    sample_sets_dict, expected_identical_diff_order, expected_partly, expected_non_overlap
):
    """
    Test to assert that different combinations of sample sets are correctly classified by the
    calculate_pairwise_overlap_types_for_sample_sets function.
    """
    result = _calculate_pairwise_overlap_types_for_sample_sets(sample_sets_dict)

    print(result)

    assert "identical elements, different order" in result
    assert "partly overlapping" in result
    assert "non-overlapping" in result

    for expected in expected_identical_diff_order:
        assert any(
            (expected[0], expected[1]) == pair or (expected[1], expected[0]) == pair
            for pair in result["identical elements, different order"]
        )

    for expected in expected_partly:
        assert any(
            (expected[0], expected[1]) == pair or (expected[1], expected[0]) == pair
            for pair in result["partly overlapping"]
        )

    for expected in expected_non_overlap:
        assert any(
            (expected[0], expected[1]) == pair or (expected[1], expected[0]) == pair
            for pair in result["non-overlapping"]
        )


@pytest.mark.parametrize(
    "files_to_download,dimensions_index,should_raise_error,expected_message_part",
    [
        # Case: identical elements, different order
        (
            ["file1.vcf.gz", "file2.vcf.gz"],
            {
                "vcf_files": [
                    {"vcf_file_s3_key": "file1.vcf.gz", "samples": ["A", "B"]},
                    {"vcf_file_s3_key": "file2.vcf.gz", "samples": ["B", "A"]},
                ]
            },
            True,
            "identical elements but different order",
        ),
        # Case: partly overlapping
        (
            ["file1.vcf.gz", "file2.vcf.gz"],
            {
                "vcf_files": [
                    {"vcf_file_s3_key": "file1.vcf.gz", "samples": ["A", "B"]},
                    {"vcf_file_s3_key": "file2.vcf.gz", "samples": ["B", "C"]},
                ]
            },
            True,
            "partly overlapping",
        ),
        # Case: no overlap (should not raise)
        (
            ["file1.vcf.gz", "file2.vcf.gz"],
            {
                "vcf_files": [
                    {"vcf_file_s3_key": "file1.vcf.gz", "samples": ["A"]},
                    {"vcf_file_s3_key": "file2.vcf.gz", "samples": ["B"]},
                ]
            },
            False,
            "No unsupported sample sets found. Proceeding with bcftools pipeline.",
        ),
    ],
)
def test_check_if_samples_can_be_combined_with_bcftools_param(
    files_to_download,
    dimensions_index,
    should_raise_error,
    expected_message_part,
    caplog,
):
    vcf_dimensions_data = dimensions_index

    if should_raise_error:
        with pytest.raises(ValueError) as excinfo:
            _check_if_samples_can_be_combined_with_bcftools(files_to_download, vcf_dimensions_data)
        assert expected_message_part in str(excinfo.value)
    else:
        with caplog.at_level("INFO"):
            _check_if_samples_can_be_combined_with_bcftools(files_to_download, vcf_dimensions_data)
        assert expected_message_part in caplog.text
