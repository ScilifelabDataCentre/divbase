from unittest.mock import MagicMock, patch

import pytest
from celery import current_app
from celery.backends.redis import RedisBackend
from kombu.connection import Connection

from divbase_lib.s3_client import create_s3_file_manager
from divbase_lib.vcf_dimension_indexing import VCFDimensionIndexManager
from divbase_worker.tasks import (
    bcftools_pipe_task,
    calculate_pairwise_overlap_types_for_sample_sets,
    check_if_samples_can_be_combined_with_bcftools,
)


@patch("divbase_worker.tasks.BcftoolsQueryManager")
@patch("divbase_lib.s3_client.S3FileManager.download_files")
@patch("divbase_lib.s3_client.S3FileManager.upload_files")
@patch(
    "divbase_worker.tasks.create_s3_file_manager",
    side_effect=lambda url=None: create_s3_file_manager(url="http://localhost:9002"),
)
@patch("divbase_worker.tasks.VCFDimensionIndexManager")
@patch("divbase_lib.queries.run_sidecar_metadata_query")
def test_bcftools_pipe_task_directly(
    mock_metadata_query,
    mock_vcf_manager,
    mock_create_s3_manager,
    mock_upload_files,
    mock_download_files,
    mock_bcftools_manager,
    bcftools_pipe_kwargs_fixture,
    tmp_path,
    sample_tsv_file,
    mock_latest_versions_of_bucket_files,
):
    """
    Run bcftools_pipe_task directly without going via a Celery worker or any other services in the job system.
    Needs many mocked services to run, so it mainly tests the task logic and flow.
    """

    mock_filenames = list(mock_latest_versions_of_bucket_files.keys())

    with patch(
        "divbase_lib.s3_client.S3FileManager.latest_version_of_all_files",
        return_value=mock_latest_versions_of_bucket_files,
    ):
        mock_vcf_manager_instance = MagicMock()
        mock_vcf_manager_instance.dimensions_info = {
            "dimensions": [
                {
                    "filename": filename,
                    "file_version_ID_in_bucket": mock_latest_versions_of_bucket_files[filename],
                    "dimensions": {
                        "variants": 10,
                        "sample_count": 5,
                        "scaffolds": sorted(["21", "24"]),
                        "sample_names": ["S1", "S2", "S3", "S4", "S5"],
                    },
                }
                for filename in mock_filenames
            ]
        }
        mock_vcf_manager_instance.get_indexed_filenames.return_value = mock_latest_versions_of_bucket_files
        mock_vcf_manager.return_value = mock_vcf_manager_instance

        dummy_vcf_files = [tmp_path / filename for filename in mock_filenames if filename.endswith(".vcf.gz")]
        mock_download_files.side_effect = [
            [sample_tsv_file],
            dummy_vcf_files,
        ]
        mock_upload_files.return_value = None

        output_file = "merged.vcf.gz"
        mock_manager_instance = MagicMock()
        mock_manager_instance.execute_pipe.return_value = output_file
        mock_bcftools_manager.return_value = mock_manager_instance

        mock_metadata_result = MagicMock()
        mock_metadata_result.unique_filenames = [
            filename for filename in mock_filenames if filename.endswith(".vcf.gz")
        ]
        mock_metadata_result.sample_and_filename_subset = [
            {"Sample_ID": "S2", "Filename": filename} for filename in mock_metadata_result.unique_filenames
        ]
        mock_metadata_result.unique_sample_ids = ["S2", "S4"]
        mock_metadata_query.return_value = mock_metadata_result

        result = bcftools_pipe_task(**bcftools_pipe_kwargs_fixture)

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
@patch("divbase_worker.tasks.BcftoolsQueryManager")
@patch("divbase_lib.s3_client.S3FileManager.download_files")
@patch("divbase_lib.queries.run_sidecar_metadata_query")
@patch("divbase_lib.s3_client.S3FileManager.upload_files")
@patch(
    "divbase_worker.tasks.create_s3_file_manager",
    side_effect=lambda url=None: create_s3_file_manager(url="http://localhost:9002"),
)
@patch("divbase_worker.tasks.VCFDimensionIndexManager")
@patch("divbase_lib.s3_client.S3FileManager.latest_version_of_all_files")
def test_bcftools_pipe_task_using_eager_mode(
    mock_latest_versions_patch,
    mock_vcf_manager,
    mock_create_s3_manager,
    mock_upload_files,
    mock_sidecar_query,
    mock_download_files,
    mock_bcftools_manager,
    bcftools_pipe_kwargs_fixture,
    sample_tsv_file,
    tmp_path,
    mock_latest_versions_of_bucket_files,
):
    """
    Test Celery task in eager mode (no worker needed).
    This is similar to running the task directly, but it uses Celery's eager mode
    which tests the task execution flow as if it were run by a worker.
    Specifically, it tests the .delay() and .get() methods of the task.
    """

    original_task_always_eager_value = current_app.conf.task_always_eager
    original_task_eager_propagates_value = current_app.conf.task_eager_propagates

    try:
        current_app.conf.update(
            task_always_eager=True,
            task_eager_propagates=True,
        )

        mock_latest_versions_patch.return_value = mock_latest_versions_of_bucket_files
        mock_filenames = list(mock_latest_versions_of_bucket_files.keys())

        mock_vcf_manager_instance = MagicMock()
        mock_vcf_manager_instance.dimensions_info = {
            "dimensions": [
                {
                    "filename": filename,
                    "file_version_ID_in_bucket": mock_latest_versions_of_bucket_files[filename],
                    "dimensions": {
                        "variants": 10,
                        "sample_count": 5,
                        "scaffolds": sorted(["21", "24"]),
                        "sample_names": ["S1", "S2", "S3", "S4", "S5"],
                    },
                }
                for filename in mock_filenames
            ]
        }
        mock_vcf_manager_instance.get_indexed_filenames.return_value = mock_latest_versions_of_bucket_files
        mock_vcf_manager.return_value = mock_vcf_manager_instance

        dummy_vcf_files = [tmp_path / filename for filename in mock_filenames if filename.endswith(".vcf.gz")]
        mock_download_files.side_effect = [
            [sample_tsv_file],
            dummy_vcf_files,
        ]
        mock_upload_files.return_value = None

        output_file = "merged.vcf.gz"
        mock_manager_instance = MagicMock()
        mock_manager_instance.execute_pipe.return_value = output_file
        mock_bcftools_manager.return_value = mock_manager_instance

        mock_metadata_result = MagicMock()
        mock_metadata_result.unique_filenames = [
            filename for filename in mock_filenames if filename.endswith(".vcf.gz")
        ]
        mock_metadata_result.sample_and_filename_subset = [
            {"Sample_ID": "S2", "Filename": filename} for filename in mock_metadata_result.unique_filenames
        ]
        mock_metadata_result.unique_sample_ids = ["S2", "S4"]
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

    This test is marked with 'integration' so it can be skipped by: pytest -m "not integration"
    Likewise, all tests marked thusly can be run with: pytest -m integration
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
