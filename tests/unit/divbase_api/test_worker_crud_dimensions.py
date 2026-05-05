"""Unit tests for worker dimensions CRUD behavior that impacts query safety."""

from unittest.mock import MagicMock

import pytest

from divbase_api.worker.crud_dimensions import get_vcf_metadata_by_project
from divbase_api.worker.tasks import _check_if_samples_can_be_combined_with_bcftools
from divbase_lib.exceptions import TaskUserError

# Shared regression guard messages
REGRESSION_GUARD_PREFIX = "Regression guard failed:"
ORDER_MUST_BE_PRESERVED_MSG = (
    f"{REGRESSION_GUARD_PREFIX} sample order must be preserved exactly from dimensions DB entries. "
    "Changing this can break bcftools concat compatibility checks downstream."
)
ORDER_NORMALIZATION_RISK_MSG = "Order normalization would hide concat-incompatible sample ordering bugs."
VALID_ORDER_UNCHANGED_MSG_TEMPLATE = (
    f"{REGRESSION_GUARD_PREFIX} expected valid sample order to remain unchanged for {{file_label}}."
)


@pytest.fixture
def mock_db():
    db = MagicMock()
    execute_result = MagicMock()
    scalars_result = MagicMock()
    execute_result.scalars.return_value = scalars_result
    db.execute.return_value = execute_result
    return db


def _make_dimensions_entry(vcf_key: str, ordered_samples: list[str]):
    """Make a mock dimensions table db entry with the given VCF key and ordered samples."""
    entry = MagicMock()
    entry.vcf_file_s3_key = vcf_key
    entry.s3_version_id = "v1"
    entry.samples = []
    for name in ordered_samples:
        sample_entry = MagicMock()
        sample_entry.sample_name = name
        entry.samples.append(sample_entry)

    scaffold = MagicMock()
    scaffold.scaffold_name = "chr1"
    entry.scaffolds = [scaffold]
    entry.variant_count = 1
    entry.sample_count = len(ordered_samples)
    entry.file_size_bytes = 1
    return entry


class TestWorkerDimensionsOrderPreservation:
    def test_regression_get_vcf_metadata_by_project_preserves_sample_order_contract(self, mock_db):
        """
        Regression test (positive outcome): dimensions retrieval must preserve sample list order from DB entries.
        Why: sample order is a hard requirement for bcftools concat compatibility checks.
        Reference: docs/development/bcftools_task_constraints.md ("Sample names must be in the same order")

        Note: scaffold order is normalized by sorting, but samples must remain in insertion order as described in the reference.
        """
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

        assert result.vcf_files[0].samples == ["S2", "S1"], ORDER_MUST_BE_PRESERVED_MSG
        assert result.vcf_files[0].scaffolds == ["chr1", "chr2"]

    def test_regression_sample_order_preservation_enforces_bcftools_concat_guard(self, mock_db):
        """
        Regression test (negative outcome): reject sample sets with identical IDs but different order.
        Rule: sample order from dimensions metadata must be preserved to ensure bcftools concat/merge can be used.
        Why: identical sample IDs in different order are concat-incompatible and must be rejected preflight.
        Reference: docs/development/bcftools_task_constraints.md ("Sample names must be in the same order")
        """

        entry_a = _make_dimensions_entry("fileA.vcf.gz", ["A", "B"])
        entry_b = _make_dimensions_entry("fileB.vcf.gz", ["B", "A"])
        mock_db.execute.return_value.scalars.return_value.all.return_value = [entry_a, entry_b]

        vcf_dimensions_data = get_vcf_metadata_by_project(db=mock_db, project_id=7)

        # Assert regression case directly: order is preserved from DB -> dimensions model
        assert vcf_dimensions_data.vcf_files[0].samples == ["A", "B"], (
            f"{REGRESSION_GUARD_PREFIX} expected sample order to remain unchanged for fileA. "
            f"{ORDER_NORMALIZATION_RISK_MSG}"
        )
        assert vcf_dimensions_data.vcf_files[1].samples == ["B", "A"], (
            f"{REGRESSION_GUARD_PREFIX} expected sample order to remain unchanged for fileB. "
            f"{ORDER_NORMALIZATION_RISK_MSG}"
        )

        with pytest.raises(TaskUserError) as excinfo:
            _check_if_samples_can_be_combined_with_bcftools(
                files_to_download=["fileA.vcf.gz", "fileB.vcf.gz"],
                vcf_dimensions_data=vcf_dimensions_data,
            )

        msg = str(excinfo.value)
        assert "Sample sets with identical elements but different order" in msg, (
            f"{REGRESSION_GUARD_PREFIX} preflight must explicitly reject identical sample IDs in different order "
            "to protect bcftools concat/merge orchestration safety."
        )
        assert "('A', 'B') vs ('B', 'A')" in msg, (
            f"{REGRESSION_GUARD_PREFIX} error message should include concrete offending sample sets for diagnosis."
        )
        assert "partly overlapping" not in msg, (
            f"{REGRESSION_GUARD_PREFIX} this case should be classified as different-order identical sets, "
            "not partly overlapping sets."
        )

    def test_regression_sample_order_preservation_allows_concat_compatible_sample_sets(self, mock_db):
        """
        Regression test (positive outcome): allow concat-compatible sample sets with identical IDs in identical order.
        Rule: sample order from dimensions metadata must be preserved for concat-compatible sample sets.
        Why: if valid identical sample order starts failing, queries get blocked even when bcftools constraints are met.
        Reference: docs/development/bcftools_task_constraints.md ("Sample names must be in the same order")
        """

        entry_a = _make_dimensions_entry("fileA.vcf.gz", ["A", "B"])
        entry_b = _make_dimensions_entry("fileB.vcf.gz", ["A", "B"])
        mock_db.execute.return_value.scalars.return_value.all.return_value = [entry_a, entry_b]

        vcf_dimensions_data = get_vcf_metadata_by_project(db=mock_db, project_id=7)

        # Assert regression case directly: valid order is preserved from DB -> dimensions model
        assert vcf_dimensions_data.vcf_files[0].samples == ["A", "B"], VALID_ORDER_UNCHANGED_MSG_TEMPLATE.format(
            file_label="fileA"
        )
        assert vcf_dimensions_data.vcf_files[1].samples == ["A", "B"], VALID_ORDER_UNCHANGED_MSG_TEMPLATE.format(
            file_label="fileB"
        )

        _check_if_samples_can_be_combined_with_bcftools(
            files_to_download=["fileA.vcf.gz", "fileB.vcf.gz"],
            vcf_dimensions_data=vcf_dimensions_data,
        )
