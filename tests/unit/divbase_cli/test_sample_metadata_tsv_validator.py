"""
Unit tests for the MetadataTSVValidator class.
"""

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


class TestDimensionMatching:
    """Test validation against project dimensions."""

    def test_samples_not_in_project(self, valid_tsv):
        """Samples in TSV but not in project should be errors."""
        project_samples = {"S1", "S2"}
        validator = MetadataTSVValidator(file_path=valid_tsv, project_samples=project_samples)
        stats, errors, warnings = validator.validate()

        assert any(
            "following samples in the TSV were not found in the DivBase project's dimensions index" in e and "S3" in e
            for e in errors
        )

    def test_samples_not_in_tsv(self, valid_tsv):
        """Samples in project but not in TSV should be warnings."""
        project_samples = {"S1", "S2", "S3", "S10", "S20"}
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
