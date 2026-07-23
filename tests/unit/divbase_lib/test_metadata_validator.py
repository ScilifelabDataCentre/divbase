"""
Test for SharedMetadataValidator: load TSV -> pandas DataFrame
-> export back to TSV and ensure that the exported file is identical to the original.

Fixtures are sourced from tests/unit/conftest.py so they are shared across all unit tests.

Fixtures intentionally EXCLUDED from identity roundtrip tests:
- sample_tsv_with_edge_cases: S2/S3 have leading/trailing whitespace in Sample_ID;
  the validator strips it, so the exported file will differ from the original.
- header_errors_tsv, sample_errors_tsv, format_errors_tsv,
  sample_tsv_with_invalid_sample_ids, sample_tsv_missing_sample_id_column,
  sample_tsv_with_duplicate_sample_ids: produce validation errors during load.
"""

from pathlib import Path

import pandas as pd
import pytest

from divbase_lib.metadata_validator import SharedMetadataValidator

FIXTURE_TSV = Path("tests/fixtures/sample_metadata.tsv")

# TODO consider if the SharedMetadataValidator and SidecarQueryManger should be updated to handle the below examples, or if they are acceptable.

# The tests in this file show that not all fixture cases used in the unit tests (tests/unit/conftest.py) will be identical after a roundtrip (TSV, load to DataFrame, export back to TSV).
# 1. numeric_multi_values_tsv: Temperature/Longitude columns mix 0 (int) with floats. Pandas writes 0 as 0.0 on export.
# 2. sample_tsv_with_numeric_data: Age column mixes 5.0 (float) with ints. Pandas writes 10 as 10.0 on export.
# 3. sample_tsv_with_edge_cases: Sample_IDs have leading/trailing whitespace that is stripped during load, so the exported file differs from the original.
ROUNDTRIP_FIXTURES = [
    "valid_tsv",
    "no_multi_values_tsv",
    "type_errors_tsv",
    "array_notation_tsv",
    "array_notation_multiple_cols_tsv",
    "sample_tsv_with_mixed_type_column",
    "sample_tsv_with_list_mixed_type_column",
]


def load_tsv(path: Path) -> pd.DataFrame:
    """Load a TSV file using SharedMetadataValidator and return the DataFrame."""
    validator = SharedMetadataValidator(file_path=path, project_samples=set(), skip_dimensions_check=True)
    result = validator.load_and_validate()
    assert result.df is not None, f"Failed to load {path}: {result.errors}"
    return result.df


def export_tsv(df: pd.DataFrame, path: Path) -> None:
    """Export DataFrame to TSV file, restoring the #Sample_ID header prefix."""
    df.copy().rename(columns={"Sample_ID": "#Sample_ID"}).to_csv(path, sep="\t", index=False)


def tsv_lines(path: Path) -> list[str]:
    """Return the non-empty lines of a TSV file for comparison."""
    return path.read_text().strip().splitlines()


class TestRoundtrip:
    """Tests that load a TSV, export it back, and verify the exported file is identical to the original."""

    @pytest.mark.parametrize("fixture_name", ROUNDTRIP_FIXTURES)
    def test_exported_tsv_identical_to_original(self, fixture_name, request, tmp_path):
        """Exported TSV content must be line-for-line identical to the original.

        Covers: basic valid data, no multi-values, negative floats, multi-value numeric cells,
        mixed-type columns (string+numeric), hyphen/range notation, array notation,
        semicolon-mixed columns, and comprehensive numeric/string data.
        """
        original: Path = request.getfixturevalue(fixture_name)
        export_path = tmp_path / "exported.tsv"
        export_tsv(load_tsv(original), export_path)
        assert tsv_lines(export_path) == tsv_lines(original)

    def test_roundtrip_no_errors(self, valid_tsv, tmp_path):
        """Exported TSV must pass validation without errors."""
        export_path = tmp_path / "exported.tsv"
        export_tsv(load_tsv(valid_tsv), export_path)
        result = SharedMetadataValidator(
            file_path=export_path, project_samples=set(), skip_dimensions_check=True
        ).load_and_validate()
        assert result.errors == [], f"Unexpected errors after roundtrip: {result.errors}"

    def test_roundtrip_fixture_tsv(self, tmp_path):
        """Real fixture TSV (tests/fixtures/sample_metadata.tsv) must be line-for-line identical after a roundtrip."""
        if not FIXTURE_TSV.exists():
            pytest.skip(f"Fixture TSV not found at {FIXTURE_TSV}")
        export_path = tmp_path / "exported.tsv"
        export_tsv(load_tsv(FIXTURE_TSV), export_path)
        assert tsv_lines(export_path) == tsv_lines(FIXTURE_TSV)
