import pandas as pd
import pytest

from divbase_tools.queries import BcftoolsQueryManager


@pytest.fixture
def bcftools_manager():
    """Return a BcftoolsQueryManager instance for testing."""
    return BcftoolsQueryManager()


@pytest.fixture
def example_sidecar_metadata_inputs_outputs():
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
def sample_tsv_file(tmp_path):
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
def tsv_with_hash_headers(tmp_path):
    """Create a TSV file with # prefix in headers."""
    data = {"#Sample_ID": ["S1", "S2"], "Population": ["Pop1", "Pop2"], "Filename": ["file1.vcf.gz", "file2.vcf.gz"]}
    df = pd.DataFrame(data)
    file_path = tmp_path / "hash_header.tsv"
    df.to_csv(file_path, sep="\t", index=False)

    return file_path
