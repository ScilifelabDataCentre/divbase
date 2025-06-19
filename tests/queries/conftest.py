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
