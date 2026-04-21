"""Unit tests for worker dimensions CRUD behavior that impacts query safety."""

from unittest.mock import MagicMock

import pytest

from divbase_api.worker.crud_dimensions import get_vcf_metadata_by_project
from divbase_api.worker.tasks import _check_if_samples_can_be_combined_with_bcftools
from divbase_lib.exceptions import TaskUserError


@pytest.fixture
def mock_db():
    db = MagicMock()
    execute_result = MagicMock()
    scalars_result = MagicMock()
    execute_result.scalars.return_value = scalars_result
    db.execute.return_value = execute_result
    return db


class TestWorkerDimensionsOrderPreservation:
    def test_get_vcf_metadata_by_project_preserves_sample_order(self, mock_db):
        """Test that worker dimensions retrieval preserves sample list order from DB entries (no sorting). Scaffolds are normalized with sorting, though."""
        entry = MagicMock()
        entry.vcf_file_s3_key = "fileA.vcf.gz"
        entry.s3_version_id = "v1"
        entry.samples = []
        for name in ["S2", "S1"]:
            sample_entry = MagicMock()
            sample_entry.sample_name = name
            entry.samples.append(sample_entry)

        entry.scaffolds = []
        for name in ["chr2", "chr1"]:
            scaffold_entry = MagicMock()
            scaffold_entry.scaffold_name = name
            entry.scaffolds.append(scaffold_entry)

        entry.variant_count = 1
        entry.sample_count = 2
        entry.file_size_bytes = 1

        entries = [entry]
        mock_db.execute.return_value.scalars.return_value.all.return_value = entries

        result = get_vcf_metadata_by_project(db=mock_db, project_id=7)

        assert result.vcf_files[0].samples == ["S2", "S1"]
        assert result.vcf_files[0].scaffolds == ["chr1", "chr2"]

    def test_order_preservation_enables_identical_elements_different_order_guard(self, mock_db):
        """
        Regression test for ensuring that sample-set order is preserved. This allows the validation logic in query manager to reject bcftools concat-incompatible files
        that have identical sample IDs but different ordering.
        """
        entry_a = MagicMock()
        entry_a.vcf_file_s3_key = "fileA.vcf.gz"
        entry_a.s3_version_id = "v1"
        entry_a.samples = []
        for name in ["A", "B"]:
            sample_entry = MagicMock()
            sample_entry.sample_name = name
            entry_a.samples.append(sample_entry)
        entry_a.scaffolds = []
        scaffold_a = MagicMock()
        scaffold_a.scaffold_name = "chr1"
        entry_a.scaffolds.append(scaffold_a)
        entry_a.variant_count = 1
        entry_a.sample_count = 2
        entry_a.file_size_bytes = 1

        entry_b = MagicMock()
        entry_b.vcf_file_s3_key = "fileB.vcf.gz"
        entry_b.s3_version_id = "v1"
        entry_b.samples = []
        for name in ["B", "A"]:
            sample_entry = MagicMock()
            sample_entry.sample_name = name
            entry_b.samples.append(sample_entry)
        entry_b.scaffolds = []
        scaffold_b = MagicMock()
        scaffold_b.scaffold_name = "chr1"
        entry_b.scaffolds.append(scaffold_b)
        entry_b.variant_count = 1
        entry_b.sample_count = 2
        entry_b.file_size_bytes = 1

        entries = [entry_a, entry_b]
        mock_db.execute.return_value.scalars.return_value.all.return_value = entries
        vcf_dimensions_data = get_vcf_metadata_by_project(db=mock_db, project_id=7)

        with pytest.raises(TaskUserError, match="identical elements but different order"):
            _check_if_samples_can_be_combined_with_bcftools(
                files_to_download=["fileA.vcf.gz", "fileB.vcf.gz"],
                vcf_dimensions_data=vcf_dimensions_data,
            )
