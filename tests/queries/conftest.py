import shutil
import tempfile
import time
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from celery import current_app

from divbase_tools.queries import BcftoolsQueryManager, SidecarQueryManager


@pytest.fixture
def bcftools_manager() -> BcftoolsQueryManager:
    """Return a BcftoolsQueryManager instance for testing."""
    return BcftoolsQueryManager()


@pytest.fixture
def create_sidecar_manager():
    def create_manager(file: Path):
        return SidecarQueryManager(file=file)

    return create_manager


@pytest.fixture
def example_sidecar_metadata_inputs_outputs() -> dict[str, Any]:
    """Return sample inputs for bcftools tests."""
    test_filenames = ["sample1.vcf.gz", "sample2.vcf.gz"]
    test_samples = [
        {"Sample_ID": "S1", "Filename": "sample1.vcf.gz"},
        {"Sample_ID": "S2", "Filename": "sample1.vcf.gz"},
        {"Sample_ID": "S3", "Filename": "sample2.vcf.gz"},
    ]
    expected_temp_files = ["temp_subset_0_0.vcf.gz", "temp_subset_0_1.vcf.gz"]
    return {
        "filenames": test_filenames,
        "sample_and_filename_subset": test_samples,
        "output_temp_files": expected_temp_files,
    }


@pytest.fixture
def sample_tsv_file(tmp_path: Path) -> Path:
    """Create a sample TSV file for testing."""
    data = {
        "Sample_ID": ["S1", "S2", "S3", "S4", "S5"],
        "Population": ["Pop1", "Pop1", "Pop2", "Pop2", "Pop3"],
        "Sex": ["M", "F", "M", "F", "M"],
        "Filename": ["file1.vcf.gz", "file1.vcf.gz", "file2.vcf.gz", "file2.vcf.gz", "file3.vcf.gz"],
    }

    df = pd.DataFrame(data)
    file_path = tmp_path / "test_samples.tsv"
    df.to_csv(file_path, sep="\t", index=False)

    return file_path


@pytest.fixture
def tsv_with_hash_headers(tmp_path: Path) -> Path:
    """Create a TSV file with # prefix in headers."""
    data = {"#Sample_ID": ["S1", "S2"], "Population": ["Pop1", "Pop2"], "Filename": ["file1.vcf.gz", "file2.vcf.gz"]}
    df = pd.DataFrame(data)
    file_path = tmp_path / "hash_header.tsv"
    df.to_csv(file_path, sep="\t", index=False)

    return file_path


@pytest.fixture
def demo_sidecar_metadata_inputs_outputs() -> dict[str, Any]:
    """Return sample inputs for bcftools tests using the VCF files and tsv query results from the initial demo."""
    test_filenames = [
        "/app/tests/fixtures/HOM_20ind_17SNPs_first_10_samples.vcf.gz",
        "/app/tests/fixtures/HOM_20ind_17SNPs_last_10_samples.vcf.gz",
    ]
    sample_and_filename_subset = [
        {"Sample_ID": "5a_HOM-I13", "Filename": test_filenames[1]},
        {"Sample_ID": "5a_HOM-I14", "Filename": test_filenames[1]},
        {"Sample_ID": "5a_HOM-I20", "Filename": test_filenames[1]},
        {"Sample_ID": "5a_HOM-I21", "Filename": test_filenames[1]},
        {"Sample_ID": "5a_HOM-I7", "Filename": test_filenames[0]},
        {"Sample_ID": "1b_HOM-G58", "Filename": test_filenames[0]},
    ]

    return {
        "filenames": test_filenames,
        "sample_and_filename_subset": sample_and_filename_subset,
    }


@pytest.fixture
def copy_fixtures_to_mock_download_from_bucket():
    """
    Fixture to simulate downloading files from a bucket by copying them from the fixtures directory.
    Designed to mock divbase_tools.services.download_files_command.
    It needs the same inputs as the download_files_command function that it is mocking, but it does not use them.
    Pytest fixtures cannot take parameters, and thus the actual function is defined and returned from inside the fixture.
    """

    def mock_download_files_command(
        file_list,
        bucket_name=None,
        all_files=None,
        download_dir=None,
        bucket_version=None,
        config_path=None,
    ):
        for file_path in file_list:
            file_name = Path(file_path).name
            destination = Path(file_name)
            if Path(file_path).exists():
                if not destination.exists():
                    shutil.copy(file_path, destination)
                    print(f"Copied {file_name} from fixtures to current directory", flush=True)
                else:
                    print(f"File {file_name} already exists in current directory, skipping copy", flush=True)
            else:
                print(f"File not found in fixtures: {file_name}")

    return mock_download_files_command


@pytest.fixture
def valid_tsv_path():
    with tempfile.NamedTemporaryFile("w", delete=False, suffix=".tsv") as tmp:
        tmp.write("Sample_ID\tFilename\tcol1\tcol2\nS1\tS1.vcf\tA\t1\nS2\tS2.vcf\tB\t2\nS3\tS3.vcf\tB\t3\n")
        tmp_path = Path(tmp.name)
    yield tmp_path
    tmp_path.unlink()


@pytest.fixture
def wait_for_celery_task_completion():
    def wait_for_task_completion(task_id: str, max_wait: int = 30) -> Any:
        """Wait for a Celery task to complete, with timeout."""
        async_result = current_app.AsyncResult(task_id)
        start_time = time.time()

        while not async_result.ready():
            if time.time() - start_time > max_wait:
                pytest.fail(f"Task timed out after {max_wait} seconds")
            time.sleep(1)

        return async_result.get()

    return wait_for_task_completion


@pytest.fixture
def concurrency_of_default_queue():
    """
    Returns a dict of {worker: concurrency} for all workers connected to the 'celery' queue.
    """
    app_stats = current_app.control.inspect().stats()
    max_concurrency_per_worker = {}
    for worker, stats in app_stats.items():
        max_concurrency_per_worker[worker] = stats.get("pool", {}).get("max-concurrency")

    queue_info = current_app.control.inspect().active_queues()
    default_queue_concurrency = {}
    for worker, queues in queue_info.items():
        for queue in queues:
            if queue["name"] == "celery":
                default_queue_concurrency[worker] = max_concurrency_per_worker[worker]
    return default_queue_concurrency


@pytest.fixture
def sample_metadata_query_kwargs_fixture():
    """Standard kwargs for testing sample_metadata_query_task."""
    return {
        "tsv_filter": "Area:West of Ireland,Northern Portugal;Sex:F",
        "metadata_tsv_name": "sample_metadata.tsv",
        "bucket_name": "query-project",
    }


@pytest.fixture
def bcftools_pipe_kwargs_fixture():
    """Standard kwargs for testing bcftools_pipe_task"""
    return {
        "tsv_filter": "Area:West of Ireland,Northern Portugal;Sex:F",
        "metadata_tsv_name": "sample_metadata.tsv",
        "command": "view -s SAMPLES; view -r 21:15000000-25000000",
        "bucket_name": "query-project",
        "user_name": "Default User",
    }
