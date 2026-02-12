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
    """Create TSV with type errors: mixed types in column and cell, hyphen in numeric, and range notation.

    Population: Has both cell-level error (1;three;5) and column-level mixed types (numeric + string)
    Test: Has column-level mixed types (all numeric values + string 'all')
    Code: String column with hyphen in one value
    Range: Contains range notation (e.g., '1-2') which should be rejected in numeric columns
    """
    tsv_content = """#Sample_ID\tPopulation\tTest\tCode\tRange
S1\t1\t2\tA100\t1-2
S2\tabc\t3\tB200\t3
S3\t1;three;5\tall\tC300\t4
S4\t3-5\t4\tD400\t5
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
    """Create a TSV file with multi-value numeric cells and negative numbers to verify they're classified as numeric."""
    tsv_content = """#Sample_ID\tScores\tValues\tTemperature\tLongitude\tLatitude\tElevation
S1\t1;2;3\t10;20\t-5.5\t-2.78305556\t51.5\t100
S2\t4;5\t30;40;50\t-10.2\t-0.12765\t52.2\t-50
S3\t6\t60\t0\t1.25\t50.8\t-100.5
S4\t7;8;9;10\t70\t15.5\t-3.5;-2.1\t49.5\t200
S5\t11\t80;90\t-20\t0\t48.2\t-25
"""
    tsv_file = tmp_path / "numeric_multi_values.tsv"
    tsv_file.write_text(tsv_content)
    return tsv_file


@pytest.fixture
def negative_numeric_columns():
    """Columns in the numeric_multi_values_tsv fixture that should be classified as numeric (including negative values)."""
    return ["Temperature", "Longitude", "Latitude", "Elevation"]


def test_valid_tsv_passes_all_checks(valid_tsv, project_samples):
    """Valid TSV should pass with no errors or warnings."""
    stats, errors, warnings = MetadataTSVValidator.validate(valid_tsv, project_samples)

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
        stats, errors, warnings = MetadataTSVValidator.validate(header_errors_tsv, project_samples)

        assert any("First column must be named '#Sample_ID'" in e for e in errors)

    def test_duplicate_column_names(self, header_errors_tsv, project_samples):
        """Duplicate column names should be detected."""
        stats, errors, warnings = MetadataTSVValidator.validate(header_errors_tsv, project_samples)

        assert any("Duplicate column names" in e and "Area" in e for e in errors)

    def test_empty_column_name(self, header_errors_tsv, project_samples):
        """Empty column names should be detected."""
        stats, errors, warnings = MetadataTSVValidator.validate(header_errors_tsv, project_samples)

        assert any("Empty column name" in e for e in errors)


class TestSampleIDValidation:
    """Test validation of Sample_ID column."""

    def test_empty_sample_id(self, sample_errors_tsv, project_samples):
        """Empty Sample_ID should be detected."""
        stats, errors, warnings = MetadataTSVValidator.validate(sample_errors_tsv, project_samples)

        assert any("Sample_ID is empty" in e for e in errors)

    def test_semicolon_in_sample_id(self, sample_errors_tsv, project_samples):
        """Sample_ID containing semicolon should be detected."""
        stats, errors, warnings = MetadataTSVValidator.validate(sample_errors_tsv, project_samples)

        assert any("contains semicolon" in e and "S3;S4" in e for e in errors)

    def test_duplicate_sample_id(self, sample_errors_tsv, project_samples):
        """Duplicate Sample_IDs should be detected."""
        stats, errors, warnings = MetadataTSVValidator.validate(sample_errors_tsv, project_samples)

        assert any("Duplicate Sample_ID" in e and "S1" in e for e in errors)


class TestFormattingValidation:
    """Test validation of TSV formatting."""

    def test_wrong_column_count(self, format_errors_tsv, project_samples):
        """Rows with wrong number of columns should be detected."""
        stats, errors, warnings = MetadataTSVValidator.validate(format_errors_tsv, project_samples)

        assert any("Expected 3 tab-separated columns" in e and "found 2" in e for e in errors)

    def test_comma_in_cell(self, format_errors_tsv, project_samples):
        """Commas in cells should be detected."""
        stats, errors, warnings = MetadataTSVValidator.validate(format_errors_tsv, project_samples)

        assert any("forbidden character ','" in e for e in errors)

    def test_whitespace_warning(self, format_errors_tsv, project_samples):
        """Leading/trailing whitespace should generate warnings."""
        stats, errors, warnings = MetadataTSVValidator.validate(format_errors_tsv, project_samples)

        assert any("leading or trailing whitespace" in w for w in warnings)


class TestTypeValidation:
    """Test validation of column types.

    Mixed types (columns with both numeric-looking and non-numeric values) are treated
    as string columns and reported as warnings, not errors.
    """

    def test_mixed_types_in_column_is_warning(self, type_errors_tsv, project_samples):
        """Columns with mixed numeric and string types should produce a warning (not error) and be classified as mixed_type."""
        stats, errors, warnings = MetadataTSVValidator.validate(type_errors_tsv, project_samples)

        assert any("mixed" in w.lower() and "Population" in w for w in warnings)
        assert not any("mixed types" in e.lower() and "Population" in e for e in errors)

    def test_mixed_types_in_cell_is_warning(self, type_errors_tsv, project_samples):
        """Cells with mixed types (e.g., '1;three;5') should produce a warning (not error)."""
        stats, errors, warnings = MetadataTSVValidator.validate(type_errors_tsv, project_samples)

        assert any("1;three;5" in w and "mixed types" in w.lower() for w in warnings)
        assert not any("1;three;5" in e and "mixed types" in e.lower() for e in errors)

    def test_hyphen_in_numeric_looking_column_is_warning(self, type_errors_tsv, project_samples):
        """Hyphens in values that look like range notation should produce a warning (not error)."""
        stats, errors, warnings = MetadataTSVValidator.validate(type_errors_tsv, project_samples)

        assert any("hyphen" in w.lower() and "3-5" in w for w in warnings)
        assert not any("hyphen" in e.lower() and "3-5" in e for e in errors)

    def test_cell_and_column_level_mixed_types_are_warnings(self, type_errors_tsv, project_samples):
        """When a column has both cell-level and column-level mixed types, both should produce warnings (not errors)."""
        stats, errors, warnings = MetadataTSVValidator.validate(type_errors_tsv, project_samples)

        assert any("1;three;5" in w and "mixed types" in w.lower() for w in warnings)
        assert any("mixed" in w.lower() and "Population" in w for w in warnings)
        assert "Population" in stats["mixed_type_columns"]
        assert "Test" in stats["mixed_type_columns"]

    def test_stats_show_mixed_type_columns(self, type_errors_tsv, project_samples):
        """
        Stats should show columns as mixed-type for informational purposes.
        The type_errors_tsv fixture used here has columns with mixed types.
        """
        stats, errors, warnings = MetadataTSVValidator.validate(type_errors_tsv, project_samples)

        assert "Population" in stats["mixed_type_columns"]
        assert "Test" in stats["mixed_type_columns"]
        assert len(stats["mixed_type_columns"]) == 3

    def test_multi_value_numeric_cells_are_numeric(self, numeric_multi_values_tsv, project_samples):
        """Multi-value numeric cells (e.g., '2;4') should be correctly classified as numeric, not string or mixed-type."""
        stats, errors, warnings = MetadataTSVValidator.validate(numeric_multi_values_tsv, project_samples)

        assert "Scores" in stats["numeric_columns"]
        assert "Values" in stats["numeric_columns"]
        assert "Scores" not in stats["string_columns"]
        assert "Values" not in stats["string_columns"]
        assert "Scores" not in stats["mixed_type_columns"]
        assert "Values" not in stats["mixed_type_columns"]
        assert not any("mixed" in w.lower() and ("Scores" in w or "Values" in w) for w in warnings)


class TestDimensionMatching:
    """Test validation against project dimensions."""

    def test_samples_not_in_project(self, valid_tsv):
        """Samples in TSV but not in project should be errors."""
        project_samples = {"S1", "S2"}
        stats, errors, warnings = MetadataTSVValidator.validate(valid_tsv, project_samples)

        assert any(
            "following samples in the TSV were not found in the DivBase project's dimensions index" in e and "S3" in e
            for e in errors
        )

    def test_samples_not_in_tsv(self, valid_tsv):
        """Samples in project but not in TSV should be warnings."""
        project_samples = {"S1", "S2", "S3", "S10", "S20"}
        stats, errors, warnings = MetadataTSVValidator.validate(valid_tsv, project_samples)

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
        stats, errors, warnings = MetadataTSVValidator.validate(valid_tsv, project_samples)

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
        stats, errors, warnings = MetadataTSVValidator.validate(no_multi_values_tsv, {"S1", "S2"})
        assert stats["has_multi_values"] is False


class TestEdgeCases:
    """Test edge cases and error conditions."""

    def test_empty_file(self, project_samples, tmp_path):
        """Empty file should be detected."""
        empty_file = tmp_path / "empty.tsv"
        empty_file.write_text("")

        stats, errors, warnings = MetadataTSVValidator.validate(empty_file, project_samples)

        assert any("File is empty" in e for e in errors)

    def test_nonexistent_file(self, project_samples):
        """Nonexistent file should be handled gracefully."""
        stats, errors, warnings = MetadataTSVValidator.validate(Path("/nonexistent/file.tsv"), project_samples)

        assert any("Failed to read file" in e for e in errors)


class TestNegativeNumbers:
    """Test that negative numbers are properly handled as numeric values."""

    def test_negative_numbers_are_numeric(self, numeric_multi_values_tsv, negative_numeric_columns):
        """Test that negative numbers are correctly classified as numeric, not flagged as errors due to hyphen check for ranges in numeric cells."""
        stats, errors, warnings = MetadataTSVValidator.validate(
            numeric_multi_values_tsv, {"S1", "S2", "S3", "S4", "S5"}
        )

        for col in negative_numeric_columns:
            assert not any("hyphen" in e.lower() and col in e for e in errors)
            assert col in stats["numeric_columns"]

        assert len(stats["mixed_type_columns"]) == 0

    def test_negative_numbers_with_semicolons(self, numeric_multi_values_tsv, negative_numeric_columns):
        """Test that negative numbers in semicolon-separated cells are handled correctly."""
        stats, errors, warnings = MetadataTSVValidator.validate(
            numeric_multi_values_tsv, {"S1", "S2", "S3", "S4", "S5"}
        )

        assert "Longitude" in negative_numeric_columns
        assert "Longitude" in stats["numeric_columns"]
        assert "Longitude" not in stats["mixed_type_columns"]
        assert not any("Longitude" in e and "mixed" in e.lower() for e in errors)

    def test_range_notation_produces_warning(self, type_errors_tsv):
        """Test that range notation like '1-2' in a mixed-type column produces a warning (column treated as string)."""
        stats, errors, warnings = MetadataTSVValidator.validate(type_errors_tsv, {"S1", "S2", "S3", "S4"})
        assert any("mixed" in w.lower() and "Range" in w for w in warnings)


class TestSemicolonColumnTypeClassification:
    """Test that the validator correctly classifies columns when semicolon-separated
    cells contain a mix of numeric and non-numeric parts."""

    @pytest.fixture
    def semicolon_mixed_tsv(self, tmp_path):
        """TSV where a column has '1;1-2' - a cell with one numeric and one non-numeric part."""
        content = "#Sample_ID\tCode\tPureNumSemicolon\n"
        content += "S1\t1;1-2\t10;20;30\n"
        content += "S2\t3\t40\n"
        content += "S3\t5\t50;60\n"
        tsv_file = tmp_path / "semicolon_mixed.tsv"
        tsv_file.write_text(content)
        return tsv_file

    def test_semicolon_cell_with_non_numeric_part_is_mixed(self, semicolon_mixed_tsv):
        """A column with cell '1;1-2' should be classified as mixed-type because '1-2' is not a number."""
        stats, errors, warnings = MetadataTSVValidator.validate(semicolon_mixed_tsv, {"S1", "S2", "S3"})

        assert "Code" in stats["mixed_type_columns"]
        assert "Code" not in stats["numeric_columns"]
        assert "Code" not in stats["string_columns"]

    def test_semicolon_cell_mixed_produces_cell_level_warning(self, semicolon_mixed_tsv):
        """A cell '1;1-2' should produce a cell-level mixed-type warning."""
        stats, errors, warnings = MetadataTSVValidator.validate(semicolon_mixed_tsv, {"S1", "S2", "S3"})

        assert any("1;1-2" in w and "mixed types" in w.lower() for w in warnings)

    def test_semicolon_cell_mixed_produces_column_level_warning(self, semicolon_mixed_tsv):
        """The column-level mixed-type warning should mention the semicolon classification rule."""
        stats, errors, warnings = MetadataTSVValidator.validate(semicolon_mixed_tsv, {"S1", "S2", "S3"})

        assert any("semicolon-separated" in w and "Code" in w for w in warnings)

    def test_semicolon_cell_mixed_is_not_error(self, semicolon_mixed_tsv):
        """Mixed types from semicolon cells should NOT produce errors."""
        stats, errors, warnings = MetadataTSVValidator.validate(semicolon_mixed_tsv, {"S1", "S2", "S3"})

        assert not any("mixed" in e.lower() for e in errors)

    def test_purely_numeric_semicolon_column_stays_numeric(self, semicolon_mixed_tsv):
        """A column with only numeric values in semicolons (e.g., '10;20;30') should be numeric."""
        stats, errors, warnings = MetadataTSVValidator.validate(semicolon_mixed_tsv, {"S1", "S2", "S3"})

        assert "PureNumSemicolon" in stats["numeric_columns"]
        assert "PureNumSemicolon" not in stats["mixed_type_columns"]
        assert "PureNumSemicolon" not in stats["string_columns"]
