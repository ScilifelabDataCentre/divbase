"""
Unit tests for SidecarQueryManager filtering
"""

import pytest

from divbase_api.services.queries import SidecarQueryManager
from divbase_lib.exceptions import SidecarInvalidFilterError


@pytest.fixture
def sample_tsv_with_numeric_data(tmp_path):
    """
    Create a temporary TSV file with numeric and string columns for testing.
    Includes semicolon-separated values in some cells. Includes both int and float
    numeric values to test that both are detected as numeric.
    """
    # Keep indentation like this to ensure that leading spaces in column 1 does not cause issues.
    tsv_content = """#Sample_ID\tPopulation\tWeight\tAge\tArea\tSingleNumber\tSingleString
S1\t1\t20.0\t5.0\tNorth\t100\tString
S2\t2;4\t25.0\t10\tEast\t200\tStrings
S3\t3\t30.0\t15\tWest;South;East\t300\tSting
S4\t4\t35.0\t20\tWest\t400\tStings
S5\t5\t40.0\t25\tNorth\t500\tThing
S6\t6\t45.0\t30\tEast\t600\tThings
S7\t1;3;5\t50.0\t35\tSouth\t700\tStrong
S8\t2\t55.0\t40\tWest\t800\tStrung
S9\t7\t62.0\t45\tNorth\t900\tStang
S10\t8\t70.0\t52\tEast\t1000\tSong
"""
    tsv_file = tmp_path / "test_metadata.tsv"
    tsv_file.write_text(tsv_content)
    return tsv_file


@pytest.fixture
def sample_tsv_with_unsupported_mixed_type_data(tmp_path):
    """
    Create a temporary TSV file with a column that containes
    mixed numeric and non-numeric values to test that this correctly raises an error.
    """
    tsv_content = """#Sample_ID\tPopulation\tWeight\tAge\tArea
S1\t1\t20.0\t5.0\tNorth
S2\t2;four;5\t25.0\t10\tEast
S3\t3\t30.0\t15\tWest;South

"""
    tsv_file = tmp_path / "test_metadata_unsupported_values.tsv"
    tsv_file.write_text(tsv_content)
    return tsv_file


class TestNumericalFilteringInequalities:
    """Test inequality operators on numerical columns."""

    def test_greater_than(self, sample_tsv_with_numeric_data):
        """Test > operator returns correct samples."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Weight:>50")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 3
        assert "S8" in sample_ids
        assert "S9" in sample_ids
        assert "S10" in sample_ids

    def test_greater_than_or_equal(self, sample_tsv_with_numeric_data):
        """Test >= operator returns correct samples."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Weight:>=50")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 4
        assert "S7" in sample_ids
        assert "S8" in sample_ids
        assert "S9" in sample_ids
        assert "S10" in sample_ids

    def test_less_than(self, sample_tsv_with_numeric_data):
        """Test < operator returns correct samples."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Age:<15")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 2
        assert "S1" in sample_ids
        assert "S2" in sample_ids

    def test_less_than_or_equal(self, sample_tsv_with_numeric_data):
        """Test <= operator returns correct samples."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Age:<=15")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 3
        assert "S1" in sample_ids
        assert "S2" in sample_ids
        assert "S3" in sample_ids

    def test_inequality_on_weight_column(self, sample_tsv_with_numeric_data):
        """Test inequality on Weight column (no semicolons, pure numeric)."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Weight:>60")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 2
        assert "S9" in sample_ids
        assert "S10" in sample_ids

    def test_inequality_on_age_column(self, sample_tsv_with_numeric_data):
        """Test inequality on Age column (no semicolons, pure numeric)."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Age:>=40")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 3
        assert "S8" in sample_ids
        assert "S9" in sample_ids
        assert "S10" in sample_ids

    def test_inequality_on_single_numeric_value_column(self, sample_tsv_with_numeric_data):
        """Test inequality on single-value numeric column (that does not have semicolon separated values)."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="SingleNumber:>600")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 4
        assert "S7" in sample_ids
        assert "S8" in sample_ids
        assert "S9" in sample_ids
        assert "S10" in sample_ids


class TestNumericalFilteringRanges:
    """Test range filtering on numerical columns."""

    def test_simple_range(self, sample_tsv_with_numeric_data):
        """Test inclusive range filtering."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Weight:30-45")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 4
        assert "S3" in sample_ids
        assert "S4" in sample_ids
        assert "S5" in sample_ids
        assert "S6" in sample_ids

    def test_range_boundaries_inclusive(self, sample_tsv_with_numeric_data):
        """Test that range boundaries are inclusive."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Age:20-30")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 3
        assert "S4" in sample_ids
        assert "S5" in sample_ids
        assert "S6" in sample_ids

    def test_range_on_weight_column(self, sample_tsv_with_numeric_data):
        """Test range filtering on Weight column (no semicolons, pure numeric)."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Weight:40-60")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 4
        assert "S5" in sample_ids
        assert "S6" in sample_ids
        assert "S7" in sample_ids
        assert "S8" in sample_ids

    def test_narrow_range(self, sample_tsv_with_numeric_data):
        """Test a narrow range returns only matching samples."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Age:20-30")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 3
        assert "S4" in sample_ids
        assert "S5" in sample_ids
        assert "S6" in sample_ids

    def test_range_on_single_numeric_value_column(self, sample_tsv_with_numeric_data):
        """Test range filtering on single-value numeric column."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="SingleNumber:350-850")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 5
        assert "S4" in sample_ids
        assert "S5" in sample_ids
        assert "S6" in sample_ids
        assert "S7" in sample_ids
        assert "S8" in sample_ids


class TestNumericalFilteringDiscreteValues:
    """Test discrete value filtering on numerical columns."""

    def test_single_discrete_value(self, sample_tsv_with_numeric_data):
        """Test filtering with a single discrete value."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Weight:50")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 1
        assert "S7" in sample_ids

    def test_multiple_discrete_values(self, sample_tsv_with_numeric_data):
        """Test filtering with multiple discrete values (OR logic)."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Weight:20,30,50")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 3
        assert "S1" in sample_ids
        assert "S3" in sample_ids
        assert "S7" in sample_ids

    def test_discrete_values_with_semicolon_separated_cells(self, sample_tsv_with_numeric_data):
        """Test discrete value filtering on Population column (string column with semicolons)."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Population:1,3,5")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 4
        assert "S1" in sample_ids
        assert "S3" in sample_ids
        assert "S5" in sample_ids
        assert "S7" in sample_ids

    def test_discrete_values_match_any_semicolon_value(self, sample_tsv_with_numeric_data):
        """Test that discrete filtering matches if ANY semicolon value matches (string matching on Population)."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Population:4")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 2
        assert "S2" in sample_ids
        assert "S4" in sample_ids

    def test_discrete_age_values(self, sample_tsv_with_numeric_data):
        """Test discrete filtering on Age column."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Age:10,25,40")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 3
        assert "S2" in sample_ids
        assert "S5" in sample_ids
        assert "S8" in sample_ids

    def test_discrete_on_single_numeric_value_column_(self, sample_tsv_with_numeric_data):
        """Test discrete value filtering on single-value numeric column."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="SingleNumber:100,600")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 2
        assert "S1" in sample_ids
        assert "S6" in sample_ids


class TestSemicolonSeparatedNumericFiltering:
    """Test that inequalities and ranges work on columns with semicolon-separated numeric values."""

    def test_inequality_on_semicolon_separated_column(self, sample_tsv_with_numeric_data):
        """Test that > operator works on Population column (semicolon-separated numbers)."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Population:>4")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 5
        assert "S5" in sample_ids
        assert "S6" in sample_ids
        assert "S7" in sample_ids
        assert "S9" in sample_ids
        assert "S10" in sample_ids

    def test_inequality_less_than_on_semicolon_separated_column(self, sample_tsv_with_numeric_data):
        """Test that < operator works on Population column."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Population:<3")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 4
        assert "S1" in sample_ids
        assert "S2" in sample_ids
        assert "S7" in sample_ids
        assert "S8" in sample_ids

    def test_range_on_semicolon_separated_column(self, sample_tsv_with_numeric_data):
        """Test that range filtering works on Population column."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Population:3-6")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 6
        assert "S2" in sample_ids
        assert "S3" in sample_ids
        assert "S4" in sample_ids
        assert "S5" in sample_ids
        assert "S6" in sample_ids
        assert "S7" in sample_ids

    def test_combined_inequality_and_discrete_on_semicolon_separated(self, sample_tsv_with_numeric_data):
        """Test combining inequality and discrete values on Population column."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Population:>6,2")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 4
        assert "S2" in sample_ids
        assert "S8" in sample_ids
        assert "S9" in sample_ids
        assert "S10" in sample_ids

    def test_range_with_semicolon_values_at_boundaries(self, sample_tsv_with_numeric_data):
        """Test that range boundaries work correctly with semicolon-separated values."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Population:1-3")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 5
        assert "S1" in sample_ids
        assert "S2" in sample_ids
        assert "S3" in sample_ids
        assert "S7" in sample_ids
        assert "S8" in sample_ids

    def test_raise_exception_on_mixed_type_column_value(self, sample_tsv_with_unsupported_mixed_type_data):
        """Test that a column with mixed numeric and non-numeric values raises an error."""
        manager = SidecarQueryManager(file=sample_tsv_with_unsupported_mixed_type_data)
        with pytest.raises(
            SidecarInvalidFilterError, match="Column 'Population' in the metadata file contains mixed types"
        ):
            manager.run_query(filter_string="Population:>2")


class TestStringColumnFiltering:
    """Test string column filtering with single and semicolon-separated values."""

    def test_single_string_value_column(self, sample_tsv_with_numeric_data):
        """Test filtering on a string column with single values (no semicolons)."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="SingleString:String")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 1
        assert "S1" in sample_ids

    def test_single_string_value_column_multiple_filters(self, sample_tsv_with_numeric_data):
        """Test filtering on a single-value string column with multiple filter values (OR logic)."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="SingleString:String,Strings")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 2
        assert "S1" in sample_ids
        assert "S2" in sample_ids

    def test_semicolon_separated_string_column(self, sample_tsv_with_numeric_data):
        """Test filtering on string column with semicolon-separated values (Area column)."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Area:West")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 3
        assert "S4" in sample_ids
        assert "S8" in sample_ids
        assert "S3" in sample_ids

    def test_semicolon_separated_string_column_any_match(self, sample_tsv_with_numeric_data):
        """Test that filtering matches if ANY semicolon-separated value matches."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Area:South")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 2
        assert "S3" in sample_ids
        assert "S7" in sample_ids
