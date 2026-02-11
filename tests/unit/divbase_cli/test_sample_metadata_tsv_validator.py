"""
Unit tests for the MetadataTSVValidator class.
"""

from pathlib import Path

import pytest

from divbase_cli.services.sample_metadata_tsv_validator import MetadataTSVValidator


@pytest.fixture
def project_samples():
    """Standard set of project samples for testing."""
    return {"S1", "S2", "S3", "S4", "S5"}


@pytest.fixture
def valid_tsv(tmp_path):
    """Create a valid TSV file that passes all validation checks and includes all project samples."""
    tsv_content = """#Sample_ID\tPopulation\tArea\tWeight
S1\t1\tNorth\t12.5
S2\t2;4\tEast\t18.8
S3\t3\tWest;South\t15.0
S4\t3;5\tSouth\t20.0
S5\t4\tNorth\t22.1
"""
    tsv_file = tmp_path / "valid.tsv"
    tsv_file.write_text(tsv_content)
    return tsv_file


@pytest.fixture
def header_errors_tsv(tmp_path):
    """Create TSV with header errors: wrong first column, duplicate columns, empty column."""
    tsv_content = """SampleID\tPopulation\tArea\tArea\t
S1\t1\tNorth\tEast\tValue
"""
    tsv_file = tmp_path / "header_errors.tsv"
    tsv_file.write_text(tsv_content)
    return tsv_file


@pytest.fixture
def sample_errors_tsv(tmp_path):
    """Create TSV with Sample_ID errors: empty, semicolons, duplicates."""
    tsv_content = """#Sample_ID\tPopulation
S1\t1
\t2
S3;S4\t3
S1\t4
"""
    tsv_file = tmp_path / "sample_errors.tsv"
    tsv_file.write_text(tsv_content)
    return tsv_file


@pytest.fixture
def format_errors_tsv(tmp_path):
    """Create TSV with formatting errors: wrong column count, commas, whitespace."""
    tsv_content = """#Sample_ID\tPopulation\tArea
S1\t1\tNorth
S2\t2,3\tEast
S3 \t 4 \t West 
S4\t5
"""
    tsv_file = tmp_path / "format_errors.tsv"
    tsv_file.write_text(tsv_content)
    return tsv_file


@pytest.fixture
def type_errors_tsv(tmp_path):
    """Create TSV with type errors: mixed types in column and cell, hyphen in numeric.

    Population: Has both cell-level error (1;three;5) and column-level mixed types (numeric + string)
    Test: Has column-level mixed types (all numeric values + string 'all')
    Code: String column with hyphen in one value
    """
    tsv_content = """#Sample_ID\tPopulation\tTest\tCode
S1\t1\t2\tA100
S2\tabc\t3\tB200
S3\t1;three;5\tall\tC300
S4\t3-5\t4\tD400
"""
    tsv_file = tmp_path / "type_errors.tsv"
    tsv_file.write_text(tsv_content)
    return tsv_file


@pytest.fixture
def no_multi_values_tsv(tmp_path):
    """Create a TSV file with no semicolon-separated values in any cell."""
    tsv_content = """#Sample_ID\tPopulation\nS1\t1\nS2\t2\n"""
    tsv_file = tmp_path / "no_multi_values.tsv"
    tsv_file.write_text(tsv_content)
    return tsv_file


@pytest.fixture
def numeric_multi_values_tsv(tmp_path):
    """Create a TSV file with multi-value numeric cells to verify they're classified as numeric."""
    tsv_content = """#Sample_ID\tScores\tValues
S1\t1;2;3\t10;20
S2\t4;5\t30;40;50
S3\t6\t60
S4\t7;8;9;10\t70
S5\t11\t80;90
"""
    tsv_file = tmp_path / "numeric_multi_values.tsv"
    tsv_file.write_text(tsv_content)
    return tsv_file


def test_valid_tsv_passes_all_checks(valid_tsv, project_samples):
    """Valid TSV should pass with no errors or warnings."""
    validator = MetadataTSVValidator(file_path=valid_tsv, project_samples=project_samples)
    stats, errors, warnings = validator.validate()

    assert len(errors) == 0
    assert len(warnings) == 0
    assert stats["total_columns"] == 4
    assert stats["user_defined_columns"] == 3
    assert stats["samples_in_tsv"] == 5
    assert stats["samples_matching_project"] == 5
    assert stats["has_multi_values"] is True
    assert "Population" in stats["numeric_columns"]
    assert "Area" in stats["string_columns"]
    assert "Weight" in stats["numeric_columns"]


class TestHeaderValidation:
    """Test validation of header row."""

    def test_wrong_first_column_name(self, header_errors_tsv, project_samples):
        """First column must be '#Sample_ID'."""
        validator = MetadataTSVValidator(file_path=header_errors_tsv, project_samples=project_samples)
        stats, errors, warnings = validator.validate()

        assert any("First column must be named '#Sample_ID'" in e for e in errors)

    def test_duplicate_column_names(self, header_errors_tsv, project_samples):
        """Duplicate column names should be detected."""
        validator = MetadataTSVValidator(file_path=header_errors_tsv, project_samples=project_samples)
        stats, errors, warnings = validator.validate()

        assert any("Duplicate column names" in e and "Area" in e for e in errors)

    def test_empty_column_name(self, header_errors_tsv, project_samples):
        """Empty column names should be detected."""
        validator = MetadataTSVValidator(file_path=header_errors_tsv, project_samples=project_samples)
        stats, errors, warnings = validator.validate()

        assert any("Empty column name" in e for e in errors)


class TestSampleIDValidation:
    """Test validation of Sample_ID column."""

    def test_empty_sample_id(self, sample_errors_tsv, project_samples):
        """Empty Sample_ID should be detected."""
        validator = MetadataTSVValidator(file_path=sample_errors_tsv, project_samples=project_samples)
        stats, errors, warnings = validator.validate()

        assert any("Sample_ID is empty" in e for e in errors)

    def test_semicolon_in_sample_id(self, sample_errors_tsv, project_samples):
        """Sample_ID containing semicolon should be detected."""
        validator = MetadataTSVValidator(file_path=sample_errors_tsv, project_samples=project_samples)
        stats, errors, warnings = validator.validate()

        assert any("contains semicolon" in e and "S3;S4" in e for e in errors)

    def test_duplicate_sample_id(self, sample_errors_tsv, project_samples):
        """Duplicate Sample_IDs should be detected."""
        validator = MetadataTSVValidator(file_path=sample_errors_tsv, project_samples=project_samples)
        stats, errors, warnings = validator.validate()

        assert any("Duplicate Sample_ID" in e and "S1" in e for e in errors)


class TestFormattingValidation:
    """Test validation of TSV formatting."""

    def test_wrong_column_count(self, format_errors_tsv, project_samples):
        """Rows with wrong number of columns should be detected."""
        validator = MetadataTSVValidator(file_path=format_errors_tsv, project_samples=project_samples)
        stats, errors, warnings = validator.validate()

        assert any("Expected 3 tab-separated columns" in e and "found 2" in e for e in errors)

    def test_comma_in_cell(self, format_errors_tsv, project_samples):
        """Commas in cells should be detected."""
        validator = MetadataTSVValidator(file_path=format_errors_tsv, project_samples=project_samples)
        stats, errors, warnings = validator.validate()

        assert any("forbidden character ','" in e for e in errors)

    def test_whitespace_warning(self, format_errors_tsv, project_samples):
        """Leading/trailing whitespace should generate warnings."""
        validator = MetadataTSVValidator(file_path=format_errors_tsv, project_samples=project_samples)
        stats, errors, warnings = validator.validate()

        assert any("leading or trailing whitespace" in w for w in warnings)


class TestTypeValidation:
    """Test validation of column types."""

    def test_mixed_types_in_column(self, type_errors_tsv, project_samples):
        """Columns with mixed numeric and string types should be detected."""
        validator = MetadataTSVValidator(file_path=type_errors_tsv, project_samples=project_samples)
        stats, errors, warnings = validator.validate()

        assert any("mixed types" in e.lower() and "Population" in e for e in errors)

    def test_mixed_types_in_cell(self, type_errors_tsv, project_samples):
        """Cells with mixed types (e.g., '1;three;5') should be detected."""
        validator = MetadataTSVValidator(file_path=type_errors_tsv, project_samples=project_samples)
        stats, errors, warnings = validator.validate()

        assert any("Cell '1;three;5' contains mixed types" in e for e in errors)

    def test_hyphen_in_numeric_column(self, type_errors_tsv, project_samples):
        """Hyphens in numeric columns should be detected."""
        validator = MetadataTSVValidator(file_path=type_errors_tsv, project_samples=project_samples)
        stats, errors, warnings = validator.validate()

        assert any("hyphen" in e.lower() and "3-5" in e for e in errors)

    def test_cell_and_column_level_mixed_types(self, type_errors_tsv, project_samples):
        """When a column has both cell-level and column-level mixed types, both errors should be reported."""
        validator = MetadataTSVValidator(file_path=type_errors_tsv, project_samples=project_samples)
        stats, errors, warnings = validator.validate()

        assert any("Cell '1;three;5' contains mixed types" in e for e in errors)
        assert any("following columns contain mixed types" in e and "Population" in e for e in errors)
        assert "Population" in stats["mixed_type_columns"]
        assert "Test" in stats["mixed_type_columns"]

    def test_stats_show_mixed_type_columns_with_cell_errors(self, type_errors_tsv, project_samples):
        """
        Stats should show columns as mixed-type even when they have cell-level errors.
        The type_errors_tsv fixture used here has two columns with mixed types.
        """
        validator = MetadataTSVValidator(file_path=type_errors_tsv, project_samples=project_samples)
        stats, errors, warnings = validator.validate()

        assert "Population" in stats["mixed_type_columns"]
        assert "Test" in stats["mixed_type_columns"]
        assert len(stats["mixed_type_columns"]) == 2

    def test_multi_value_numeric_cells_are_numeric(self, numeric_multi_values_tsv, project_samples):
        """Multi-value numeric cells (e.g., '2;4') should be correctly classified as numeric, not string or mixed-type."""
        validator = MetadataTSVValidator(file_path=numeric_multi_values_tsv, project_samples=project_samples)
        stats, errors, warnings = validator.validate()

        assert "Scores" in stats["numeric_columns"]
        assert "Values" in stats["numeric_columns"]
        assert "Scores" not in stats["string_columns"]
        assert "Values" not in stats["string_columns"]
        assert "Scores" not in stats["mixed_type_columns"]
        assert "Values" not in stats["mixed_type_columns"]
        assert not any("mixed types" in e.lower() and ("Scores" in e or "Values" in e) for e in errors)


class TestDimensionMatching:
    """Test validation against project dimensions."""

    def test_samples_not_in_project(self, valid_tsv):
        """Samples in TSV but not in project should be errors."""
        project_samples = {"S1", "S2"}  # Only S1 and S2 exist in project dimensions
        validator = MetadataTSVValidator(file_path=valid_tsv, project_samples=project_samples)
        stats, errors, warnings = validator.validate()

        assert any(
            "following samples in the TSV were not found in the DivBase project's dimensions index" in e and "S3" in e
            for e in errors
        )

    def test_samples_not_in_tsv(self, valid_tsv):
        """Samples in project but not in TSV should be warnings."""
        project_samples = {"S1", "S2", "S3", "S10", "S20"}  # S10 and S20 not in TSV
        validator = MetadataTSVValidator(file_path=valid_tsv, project_samples=project_samples)
        stats, errors, warnings = validator.validate()

        assert any(
            "following samples in the DivBase project's dimensions index were not found in the TSV" in w and "S10" in w
            for w in warnings
        )
        assert any(
            "following samples in the DivBase project's dimensions index were not found in the TSV" in w and "S20" in w
            for w in warnings
        )


class TestStatistics:
    """Test statistics collection."""

    def test_statistics_collection(self, valid_tsv, project_samples):
        """Verify statistics are correctly collected."""
        validator = MetadataTSVValidator(file_path=valid_tsv, project_samples=project_samples)
        stats, errors, warnings = validator.validate()

        assert stats["total_columns"] == 4
        assert stats["user_defined_columns"] == 3
        assert stats["samples_in_tsv"] == 5
        assert stats["samples_matching_project"] == 5
        assert stats["total_project_samples"] == 5
        assert len(stats["numeric_columns"]) == 2
        assert len(stats["string_columns"]) == 1
        assert len(stats["mixed_type_columns"]) == 0
        assert stats["has_multi_values"] is True

    def test_no_multi_values_detected(self, no_multi_values_tsv):
        """Test detection when no semicolon-separated values present."""
        validator = MetadataTSVValidator(file_path=no_multi_values_tsv, project_samples={"S1", "S2"})
        stats, errors, warnings = validator.validate()
        assert stats["has_multi_values"] is False


class TestEdgeCases:
    """Test edge cases and error conditions."""

    def test_empty_file(self, project_samples, tmp_path):
        """Empty file should be detected."""
        empty_file = tmp_path / "empty.tsv"
        empty_file.write_text("")

        validator = MetadataTSVValidator(file_path=empty_file, project_samples=project_samples)
        stats, errors, warnings = validator.validate()

        assert any("File is empty" in e for e in errors)

    def test_nonexistent_file(self, project_samples):
        """Nonexistent file should be handled gracefully."""
        validator = MetadataTSVValidator(file_path=Path("/nonexistent/file.tsv"), project_samples=project_samples)
        stats, errors, warnings = validator.validate()

        assert any("Failed to read file" in e for e in errors)
