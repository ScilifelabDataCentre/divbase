from unittest.mock import MagicMock, patch

import pytest
from celery import current_app
from celery.backends.redis import RedisBackend
from kombu.connection import Connection

from divbase_api.worker.tasks import (
    bcftools_pipe_task,
    calculate_pairwise_overlap_types_for_sample_sets,
    check_if_samples_can_be_combined_with_bcftools,
)
from divbase_lib.vcf_dimension_indexing import VCFDimensionIndexManager


@pytest.mark.integration
def test_bcftools_pipe_task_with_real_worker(
    wait_for_celery_task_completion, bcftools_pipe_kwargs_fixture, run_update_dimensions
):
    """
    Integration test in which bcftools_pipe_task is run with a real Celery worker.
    Runs locally using the docker compose testing stack defined and performs a real sidecar and bcftools
    query by loading VCF files from the tests/fixtures dir. It was designed for having RabbitMQ as the broker, Redis as the backend,
    and a custom Celery worker image that has bcftools installed.
    (this test does not download fixture from bucket, since that is handled by the CLI layer)

    Does not assert file download and uploads, since that is handled by tests in tests/cli_commands/test_query_cli.py.
    (but the current version of the task does I/O with files from the bucket, so it is tested indirectly)

    """

    bucket_name = bcftools_pipe_kwargs_fixture["bucket_name"]
    run_update_dimensions(bucket_name=bucket_name)

    broker_url = current_app.conf.broker_url
    with Connection(broker_url) as conn:
        conn.ensure_connection(max_retries=1)

    if isinstance(current_app.backend, RedisBackend):
        current_app.backend.client.ping()

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
    result = calculate_pairwise_overlap_types_for_sample_sets(sample_sets_dict)

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
                "dimensions": [
                    {"filename": "file1.vcf.gz", "dimensions": {"sample_names": ["A", "B"]}},
                    {"filename": "file2.vcf.gz", "dimensions": {"sample_names": ["B", "A"]}},
                ]
            },
            True,
            "identical elements but different order",
        ),
        # Case: partly overlapping
        (
            ["file1.vcf.gz", "file2.vcf.gz"],
            {
                "dimensions": [
                    {"filename": "file1.vcf.gz", "dimensions": {"sample_names": ["A", "B"]}},
                    {"filename": "file2.vcf.gz", "dimensions": {"sample_names": ["B", "C"]}},
                ]
            },
            True,
            "partly overlapping",
        ),
        # Case: no overlap (should not raise)
        (
            ["file1.vcf.gz", "file2.vcf.gz"],
            {
                "dimensions": [
                    {"filename": "file1.vcf.gz", "dimensions": {"sample_names": ["A"]}},
                    {"filename": "file2.vcf.gz", "dimensions": {"sample_names": ["B"]}},
                ]
            },
            False,
            "No unsupported sample sets found. Proceeding with bcftools pipeline.",
        ),
    ],
)
@patch("divbase_lib.vcf_dimension_indexing.VCFDimensionIndexManager._get_bucket_dimensions_file")
def test_check_if_samples_can_be_combined_with_bcftools_param(
    mock_get_dimensions_file,
    files_to_download,
    dimensions_index,
    should_raise_error,
    expected_message_part,
    caplog,
):
    mock_get_dimensions_file.return_value = dimensions_index

    s3_file_manager = MagicMock()
    vcf_dimensions_manager = VCFDimensionIndexManager("dummy-bucket", s3_file_manager)

    if should_raise_error:
        with pytest.raises(ValueError) as excinfo:
            check_if_samples_can_be_combined_with_bcftools(files_to_download, vcf_dimensions_manager)
        assert expected_message_part in str(excinfo.value)
    else:
        with caplog.at_level("INFO"):
            check_if_samples_can_be_combined_with_bcftools(files_to_download, vcf_dimensions_manager)
        assert expected_message_part in caplog.text
