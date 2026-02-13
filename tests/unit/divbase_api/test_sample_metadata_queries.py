"""
Unit tests for SidecarQueryManager filtering
"""

import pandas as pd
import pytest

from divbase_api.services.queries import SidecarQueryManager
from divbase_lib.exceptions import (
    SidecarColumnNotFoundError,
    SidecarMetadataFormatError,
    SidecarSampleIDError,
)


@pytest.fixture
def sample_tsv_with_numeric_data(tmp_path):
    """
    Create a temporary TSV file with numeric and string columns for testing.
    Includes semicolon-separated values in some cells. Includes both int and float
    numeric values to test that both are detected as numeric. Also includes negative
    numbers to verify they are properly handled as numeric values.
    """
    # Keep indentation like this to ensure that leading spaces in column 1 does not cause issues.
    tsv_content = """#Sample_ID\tPopulation\tWeight\tAge\tArea\tSingleNumber\tSingleString\tTemperature\tLongitude\tLatitude\tElevation
S1\t1\t20.2\t5.0\tNorth\t100\tString\t-5.5\t-2.78305556\t51.5\t100
S2\t2;4\t25.0\t10\tEast\t200\tStrings\t-10.2\t-0.12765\t52.2\t-50
S3\t3\t30.8\t15\tWest;South;East\t300\tSting\t0\t1.25\t50.8\t-100.5
S4\t4\t35.1\t20\tWest\t400\tStings\t15.5\t-3.5;-2.1\t49.5\t200
S5\t5\t40.0\t25\tNorth\t500\tThing\t-20\t0\t48.2\t-25
S6\t6\t45.4\t30\tEast\t600\tThings\t10\t2.5\t53.1\t150
S7\t1;3;5\t50.9\t35\tSouth\t700\tStrong\t5\t-1.5\t52.8\t50
S8\t2\t55.2\t40\tWest\t800\tStrung\t20\t3.0\t51.0\t75
S9\t7\t62.6\t45\tNorth\t900\tStang\t-15\t-2.0\t54.5\t-10
S10\t8\t70.7\t52\tEast\t1000\tSong\t25\t1.5\t50.5\t200
"""
    tsv_file = tmp_path / "test_metadata.tsv"
    tsv_file.write_text(tsv_content)
    return tsv_file


@pytest.fixture
def sample_tsv_with_edge_cases(tmp_path):
    """
    Create a temporary TSV file to test edge cases:
    1. "string;string;string" - OK (pure strings)
    2. "1;two;5" - Mixed numeric and non-numeric: treated as string column (not an error)
    3. String values containing numbers like "1string" - OK (inferred as string)
    4. Unicode strings with diacritics - OK

    Note that S2 and S3 have leading/trailing whitespace in the Sample_ID and the code should handle that by stripping whitespace.
    """
    tsv_content = """#Sample_ID\tPureStrings\tMixedTypes\tSingleString\tSingleNumber\tUnicodeStrings\tStringWithHyphen\tNumericalWithHyphen
S1\tNorth;South;East\t1;two;5\tWest\t100\tStockholm;Göteborg\tNorth-East\t1-2
S2 \tWest;East;North\t2;three;6\tNorth\t200\tMalmö;Uppsala\tSouth-West\t2-3
 S3\tSouth\t3\tEast\t300\tKöpenhamn;København\tNorth-North-West\t3-4
S4\t1string\tstring4\tString5\t400\tHumlebæk\tEast-South-East\t4-5
"""
    tsv_file = tmp_path / "test_metadata_edge_cases.tsv"
    tsv_file.write_text(tsv_content)
    return tsv_file


@pytest.fixture
def sample_tsv_with_invalid_sample_ids(tmp_path):
    """
    Create a temporary TSV file to test Sample_ID validation:
    Has empty and duplicate Sample_IDs that both should raise error during load
    """
    tsv_content = """#Sample_ID\tPopulation\tWeight
S1\t1\t20.2
\t2\t25.0
S3\t3\t30.8
S3\t4\t35.1
"""
    tsv_file = tmp_path / "test_metadata_invalid_sample_ids.tsv"
    tsv_file.write_text(tsv_content)
    return tsv_file


@pytest.fixture
def sample_tsv_missing_sample_id_column(tmp_path):
    """
    Create a temporary TSV file that omits the Sample_ID column.
    Should trigger SidecarColumnNotFoundError during file load.
    """
    tsv_content = """Population\tWeight\tAge\tArea
1\t20.2\t5.0\tNorth
2\t25.0\t10\tEast
3\t30.8\t15\tWest
"""
    tsv_file = tmp_path / "test_metadata_missing_sample_id.tsv"
    tsv_file.write_text(tsv_content)
    return tsv_file


@pytest.fixture
def sample_tsv_with_duplicate_sample_ids(tmp_path):
    """
    Create a temporary TSV file to test duplicate Sample_IDs (should raise error during load).
    """
    tsv_content = """#Sample_ID\tPopulation\tWeight
S1\t1\t20.2
S2\t2\t25.0
S3\t3\t30.8
S3\t4\t35.1
S4\t5\t40.0
"""
    tsv_file = tmp_path / "test_duplicate_sample_ids.tsv"
    tsv_file.write_text(tsv_content)
    return tsv_file


@pytest.fixture
def sample_tsv_with_mixed_type_column(tmp_path):
    """
    Create a TSV with a column that has mixed numeric-looking and non-numeric values,
    similar to Population_code with values like "8", "1a", "5a".
    For testing query warnings for mixed-type columns.
    """
    tsv_content = """#Sample_ID\tPopulation_code\tArea\tWeight
S1\t8\tNorth\t12.5
S2\t1a\tEast\t18.8
S3\t5a\tWest\t15.0
S4\t1b\tSouth\t20.0
S5\t4\tNorth\t22.1
"""
    tsv_file = tmp_path / "test_mixed_type.tsv"
    tsv_file.write_text(tsv_content)
    return tsv_file


@pytest.fixture
def sample_tsv_with_semicolon_mixed_type_column(tmp_path):
    """
    Create a TSV where a column has semicolon-separated values with a non-numeric part,
    e.g. '1;1-2'. This makes the column a string column because '1-2' is not a number.
    Tests that semicolon-split values are individually checked for numeric parsing.
    """
    tsv_content = """#Sample_ID\tCode\tPureNumericSemicolon\tWeight
S1\t1;1-2\t10;20;30\t12.5
S2\t3\t40\t18.8
S3\t5\t50;60\t15.0
S4\t7\t70;80;90\t20.0
"""
    tsv_file = tmp_path / "test_semicolon_mixed.tsv"
    tsv_file.write_text(tsv_content)
    return tsv_file


class TestNumericalFilteringInequalities:
    """Test inequality operators on numerical columns."""

    def test_greater_than(self, sample_tsv_with_numeric_data):
        """Test > operator returns correct samples."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Weight:>50")

        sample_ids = result.get_unique_values("Sample_ID")
        # Weight > 50: S7 (50.9), S8 (55.2), S9 (62.6), S10 (70.7)
        assert len(sample_ids) == 4
        assert "S7" in sample_ids
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

    def test_not_operator_with_inequality(self, sample_tsv_with_numeric_data):
        """Test NOT operator (!) with inequality: Population:<4,!2 should return 1 and 3 but not 2."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Population:<4,!2")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 3
        assert "S1" in sample_ids
        assert "S3" in sample_ids
        assert "S7" in sample_ids
        assert "S2" not in sample_ids
        assert "S8" not in sample_ids

    def test_not_operator_standalone(self, sample_tsv_with_numeric_data):
        """Test NOT operator (!) standalone: Population:!2 should return all except 2."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Population:!2")

        sample_ids = result.get_unique_values("Sample_ID")
        assert "S1" in sample_ids
        assert "S3" in sample_ids
        assert "S4" in sample_ids
        assert "S5" in sample_ids
        assert "S6" in sample_ids
        assert "S7" in sample_ids
        assert "S9" in sample_ids
        assert "S10" in sample_ids
        assert "S2" not in sample_ids
        assert "S8" not in sample_ids


class TestNumericalFilteringRanges:
    """Test range filtering on numerical columns."""

    def test_simple_range(self, sample_tsv_with_numeric_data):
        """Test inclusive range filtering."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Weight:30-45")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 3
        assert "S3" in sample_ids
        assert "S4" in sample_ids
        assert "S5" in sample_ids

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

    def test_not_operator_with_range(self, sample_tsv_with_numeric_data):
        """Test NOT operator (!) with range: Age:!20-30 should exclude values in range 20-30."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Age:!20-30")

        sample_ids = result.get_unique_values("Sample_ID")
        assert "S1" in sample_ids
        assert "S2" in sample_ids
        assert "S3" in sample_ids
        assert "S7" in sample_ids
        assert "S8" in sample_ids
        assert "S9" in sample_ids
        assert "S10" in sample_ids
        assert "S4" not in sample_ids
        assert "S5" not in sample_ids
        assert "S6" not in sample_ids


class TestNumericalFilteringDiscreteValues:
    """Test discrete value filtering on numerical columns."""

    def test_single_discrete_value(self, sample_tsv_with_numeric_data):
        """Test filtering with a single discrete value."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Weight:50.9")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 1
        assert "S7" in sample_ids

    def test_multiple_discrete_values(self, sample_tsv_with_numeric_data):
        """Test filtering with multiple discrete values (OR logic)."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Weight:20.2,30.8,50.9")

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

    def test_not_operator_with_discrete_values(self, sample_tsv_with_numeric_data):
        """Test NOT operator (!) with discrete values: Population:1,3,!2 should return 1 and 3 but not 2."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Population:1,3,!2")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 3
        assert "S1" in sample_ids
        assert "S3" in sample_ids
        assert "S7" in sample_ids
        assert "S2" not in sample_ids
        assert "S8" not in sample_ids

    def test_not_operator_multiple_negations(self, sample_tsv_with_numeric_data):
        """Test NOT operator (!) with multiple negations: Weight:>30,!50.9,!55.2 should exclude 50.9 and 55.2."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Weight:>30,!50.9,!55.2")

        sample_ids = result.get_unique_values("Sample_ID")
        assert "S4" in sample_ids
        assert "S5" in sample_ids
        assert "S6" in sample_ids
        assert "S9" in sample_ids
        assert "S10" in sample_ids
        assert "S7" not in sample_ids
        assert "S8" not in sample_ids
        assert "S8" not in sample_ids


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

    def test_not_operator_with_semicolon_separated(self, sample_tsv_with_numeric_data):
        """Test NOT operator (!) with semicolon-separated values: Population:>3,!5 should exclude 5."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Population:>3,!5")

        sample_ids = result.get_unique_values("Sample_ID")
        assert "S2" in sample_ids
        assert "S4" in sample_ids
        assert "S6" in sample_ids
        assert "S9" in sample_ids
        assert "S10" in sample_ids
        assert "S5" not in sample_ids
        assert "S7" not in sample_ids


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

    def test_not_operator_with_string_values(self, sample_tsv_with_numeric_data):
        """Test NOT operator (!) with string values: Area:!North should exclude North."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Area:!North")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 7
        assert "S2" in sample_ids
        assert "S3" in sample_ids
        assert "S4" in sample_ids
        assert "S6" in sample_ids
        assert "S7" in sample_ids
        assert "S8" in sample_ids
        assert "S10" in sample_ids
        assert "S1" not in sample_ids
        assert "S5" not in sample_ids
        assert "S9" not in sample_ids

    def test_not_operator_with_string_positive_and_negative(self, sample_tsv_with_numeric_data):
        """Test NOT operator (!) combined with positive values: Area:East,West,!South."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Area:East,West,!South")

        sample_ids = result.get_unique_values("Sample_ID")
        assert "S2" in sample_ids
        assert "S4" in sample_ids
        assert "S6" in sample_ids
        assert "S8" in sample_ids
        assert "S10" in sample_ids
        assert "S3" not in sample_ids
        assert "S7" not in sample_ids


class TestEdgeCases:
    """Edge case tests for SidecarQueryManager filtering."""

    def test_mixed_types_treated_as_string(self, sample_tsv_with_edge_cases):
        """Test that a column with mixed numeric and non-numeric values is treated as a string column.
        The MixedTypes column has values like '1;two;5', '2;three;6', '3', 'string4'.
        When treated as string, filtering for '1' should match cells containing '1' as a semicolon-separated value."""
        manager = SidecarQueryManager(file=sample_tsv_with_edge_cases)
        result = manager.run_query(filter_string="MixedTypes:1")
        sample_ids = result.get_unique_values("Sample_ID")
        assert "S1" in sample_ids

    def test_string_with_numbers_in_value(self, sample_tsv_with_edge_cases):
        """Test that values like '1string' are correctly inferred as strings."""
        manager = SidecarQueryManager(file=sample_tsv_with_edge_cases)
        result = manager.run_query(filter_string="PureStrings:1string")
        sample_ids = result.get_unique_values("Sample_ID")
        assert "S4" in sample_ids

    def test_unicode_string_filtering(self, sample_tsv_with_edge_cases):
        """Test that filtering for Unicode values like 'Göteborg' and 'Malmö' works and returns correct samples."""
        manager = SidecarQueryManager(file=sample_tsv_with_edge_cases)
        result = manager.run_query(filter_string="UnicodeStrings:Göteborg")
        sample_ids = result.get_unique_values("Sample_ID")
        assert "S1" in sample_ids

    def test_multi_column_single_string_and_single_number(self, sample_tsv_with_edge_cases):
        """Test that filtering on two valid single-value columns (SingleString and SingleNumber) will pass."""
        manager = SidecarQueryManager(file=sample_tsv_with_edge_cases)
        result = manager.run_query(filter_string="SingleString:String5;SingleNumber:400")
        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 1
        assert "S4" in sample_ids

    def test_multi_column_single_string_and_mixed_types_treated_as_string(self, sample_tsv_with_edge_cases):
        """Test that filtering on SingleString and MixedTypes (treated as string) works with string matching."""
        manager = SidecarQueryManager(file=sample_tsv_with_edge_cases)
        result = manager.run_query(filter_string="SingleString:String5;MixedTypes:string4")
        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 1
        assert "S4" in sample_ids

    def test_multi_column_single_number_and_mixed_types_treated_as_string(self, sample_tsv_with_edge_cases):
        """Test that filtering on SingleNumber (numeric) and MixedTypes (treated as string) works correctly."""
        manager = SidecarQueryManager(file=sample_tsv_with_edge_cases)
        result = manager.run_query(filter_string="SingleNumber:400;MixedTypes:string4")
        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 1
        assert "S4" in sample_ids

    def test_multi_column_single_string_and_pure_strings(self, sample_tsv_with_edge_cases):
        """Test that filtering on SingleString and PureStrings (both valid string columns) will pass."""
        manager = SidecarQueryManager(file=sample_tsv_with_edge_cases)
        result = manager.run_query(filter_string="SingleString:String5;PureStrings:1string")
        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 1
        assert "S4" in sample_ids

    def test_multi_column_with_unicode(self, sample_tsv_with_edge_cases):
        """Test that multi-column filtering works with unicode strings."""
        manager = SidecarQueryManager(file=sample_tsv_with_edge_cases)
        result = manager.run_query(filter_string="UnicodeStrings:København;SingleString:East")
        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 1
        assert "S3" in sample_ids

    def test_hyphens_allowed_in_string_values(self, sample_tsv_with_edge_cases):
        """Test that hyphens are allowed in string values and can be filtered correctly."""
        manager = SidecarQueryManager(file=sample_tsv_with_edge_cases)
        result = manager.run_query(filter_string="StringWithHyphen:South-West")
        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 1
        assert "S2" in sample_ids

    def test_hyphens_in_column_treated_as_string(self, sample_tsv_with_edge_cases):
        """Test that a column with hyphenated values like '1-2', '2-3' is treated as a string column.
        The values are not parseable as floats, so the column is inferred as string. Querying for
        the exact string value '2-3' should return the matching row."""
        manager = SidecarQueryManager(file=sample_tsv_with_edge_cases)
        result = manager.run_query(filter_string="NumericalWithHyphen:2-3")
        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 1
        assert "S2" in sample_ids

    def test_not_operator_edge_case_with_unicode(self, sample_tsv_with_edge_cases):
        """Test NOT operator (!) with unicode string values."""
        manager = SidecarQueryManager(file=sample_tsv_with_edge_cases)
        result = manager.run_query(filter_string="UnicodeStrings:!København")
        sample_ids = result.get_unique_values("Sample_ID")
        assert "S1" in sample_ids
        assert "S2" in sample_ids
        assert "S4" in sample_ids
        assert "S3" not in sample_ids

    def test_not_operator_only_negations(self, sample_tsv_with_edge_cases):
        """Test NOT operator (!) with only negations (no positive values): PureStrings:!North should return all except North."""
        manager = SidecarQueryManager(file=sample_tsv_with_edge_cases)
        result = manager.run_query(filter_string="PureStrings:!North")
        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 2
        assert "S3" in sample_ids
        assert "S4" in sample_ids
        assert "S1" not in sample_ids
        assert "S2" not in sample_ids


class TestSampleIDValidation:
    """Test Sample_ID validation during file loading."""

    def test_empty_sample_id_raises_error(self, sample_tsv_with_invalid_sample_ids):
        """Test that empty Sample_ID values raise SidecarSampleIDError directly during file load."""
        with pytest.raises(SidecarSampleIDError) as excinfo:
            SidecarQueryManager(file=sample_tsv_with_invalid_sample_ids)
        assert "Sample_ID column contains empty or missing values" in str(excinfo.value)

    def test_duplicate_sample_id_raises_error(self, sample_tsv_with_duplicate_sample_ids):
        """Test that duplicate Sample_ID values raise SidecarSampleIDError directly during file load."""
        with pytest.raises(SidecarSampleIDError) as excinfo:
            SidecarQueryManager(file=sample_tsv_with_duplicate_sample_ids)
        assert "Duplicate Sample_IDs found" in str(excinfo.value)
        assert "S3" in str(excinfo.value)

    def test_missing_sample_id_column_raises_error(self, sample_tsv_missing_sample_id_column):
        """Test that missing Sample_ID column raises SidecarColumnNotFoundError during file load."""
        with pytest.raises(SidecarColumnNotFoundError) as excinfo:
            SidecarQueryManager(file=sample_tsv_missing_sample_id_column)
        assert "The 'Sample_ID' column is required in the metadata file." in str(excinfo.value)


class TestNegativeNumbers:
    """Test that negative numbers are properly handled as numeric values."""

    def test_negative_numbers_in_single_value_column(self, sample_tsv_with_numeric_data):
        """Test that negative numbers in single-value columns are treated as numeric and can be filtered with inequalities."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Temperature:<0")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 4
        assert "S1" in sample_ids
        assert "S2" in sample_ids
        assert "S5" in sample_ids
        assert "S9" in sample_ids

    def test_negative_numbers_discrete_values(self, sample_tsv_with_numeric_data):
        """Test that negative numbers can be used as discrete filter values."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Temperature:-5.5,-20")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 2
        assert "S1" in sample_ids
        assert "S5" in sample_ids

    def test_negative_numbers_greater_than_inequality(self, sample_tsv_with_numeric_data):
        """Test greater than with negative numbers."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Temperature:<-5")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 4
        assert "S1" in sample_ids
        assert "S2" in sample_ids
        assert "S5" in sample_ids
        assert "S9" in sample_ids

    def test_negative_numbers_in_semicolon_cells(self, sample_tsv_with_numeric_data):
        """Test that negative numbers in semicolon-separated cells work correctly."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Longitude:-3.5")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 1
        assert "S4" in sample_ids


class TestQueryWarnings:
    """Test that the query engine produces helpful warnings when users
    filter on mixed-type or string columns with numeric syntax."""

    def test_mixed_type_column_warns_on_inequality_filter(self, sample_tsv_with_mixed_type_column):
        """Test that filtering with inequality syntax on a mixed-type column should produce a warning
        that mentions both mixed types and the comparison operators."""
        manager = SidecarQueryManager(file=sample_tsv_with_mixed_type_column)
        result = manager.run_query(filter_string="Population_code:>5")

        assert any(
            "mixed types" in w.lower() and "comparison operators" in w.lower() and "Population_code" in w
            for w in result.warnings
        )

    def test_mixed_type_column_range_syntax_does_string_match(self, sample_tsv_with_mixed_type_column):
        """Test that range-like patterns (e.g., '1-5') on a mixed-type column are treated as literal string
        matches, not numeric ranges. A general mixed-type warning is expected, but not a
        'numeric operations won't work' warning since hyphenated values are common in strings."""
        manager = SidecarQueryManager(file=sample_tsv_with_mixed_type_column)
        result = manager.run_query(filter_string="Population_code:1-5")

        assert any("mixed types" in w.lower() and "Population_code" in w for w in result.warnings)
        assert not any("will not work" in w for w in result.warnings)

    def test_mixed_type_column_no_warning_on_string_filter(self, sample_tsv_with_mixed_type_column):
        """Test that filtering with plain string values on a mixed-type column should produce a
        general mixed-type info warning but not a numeric-syntax warning."""
        manager = SidecarQueryManager(file=sample_tsv_with_mixed_type_column)
        result = manager.run_query(filter_string="Population_code:8,1a")

        sample_ids = result.get_unique_values("Sample_ID")
        assert "S1" in sample_ids
        assert "S2" in sample_ids
        assert any("mixed types" in w.lower() and "Population_code" in w for w in result.warnings)
        assert not any("will not work" in w for w in result.warnings)

    def test_pure_string_column_warns_on_numeric_inequality_filter(self, sample_tsv_with_mixed_type_column):
        """Test that filtering with numeric inequality syntax (>5) on a pure string column (Area) should result in a warning."""
        manager = SidecarQueryManager(file=sample_tsv_with_mixed_type_column)
        result = manager.run_query(filter_string="Area:>5")

        assert any(
            "string column" in w.lower() and "Area" in w and "comparison operators" in w.lower()
            for w in result.warnings
        )
        assert not any("mixed types" in w.lower() and "Area" in w for w in result.warnings)

    def test_pure_string_column_no_warning_on_normal_filter(self, sample_tsv_with_mixed_type_column):
        """Test that filtering with plain string values on a pure string column should produce no warnings."""
        manager = SidecarQueryManager(file=sample_tsv_with_mixed_type_column)
        result = manager.run_query(filter_string="Area:North,East")

        sample_ids = result.get_unique_values("Sample_ID")
        assert "S1" in sample_ids
        assert "S2" in sample_ids
        assert len(result.warnings) == 0

    def test_numeric_column_no_false_warning(self, sample_tsv_with_mixed_type_column):
        """Test that filtering with numeric syntax on a numeric column should NOT produce a warning."""
        manager = SidecarQueryManager(file=sample_tsv_with_mixed_type_column)
        result = manager.run_query(filter_string="Weight:>15")

        sample_ids = result.get_unique_values("Sample_ID")
        assert "S2" in sample_ids
        assert "S4" in sample_ids
        assert "S5" in sample_ids
        assert not any("string column" in w.lower() for w in result.warnings)
        assert not any("comparison operators" in w.lower() for w in result.warnings)

    def test_warning_mentions_semicolon_rule(self, sample_tsv_with_mixed_type_column):
        """Test that query warnings explain the semicolon classification rule."""
        manager = SidecarQueryManager(file=sample_tsv_with_mixed_type_column)
        result = manager.run_query(filter_string="Population_code:>5")

        assert any("semicolon-separated" in w for w in result.warnings)

    @pytest.mark.parametrize(
        "column,filter_string,expected_warning,expected_sample_ids",
        [
            ("Area", ">North", ["comparison operators", "Area", "only work on numeric"], []),
            ("Area", "<East", ["comparison operators", "Area"], None),
            ("Area", ">=North", ["comparison operators", ">=North"], None),
            ("Population_code", ">8", ["mixed types", "comparison operators"], None),
        ],
    )
    def test_comparison_operator_parametrized(
        self, sample_tsv_with_mixed_type_column, column, filter_string, expected_warning, expected_sample_ids
    ):
        """
        Test for comparison operator warnings and result length on string and mixed-type columns.
        """
        manager = SidecarQueryManager(file=sample_tsv_with_mixed_type_column)
        result = manager.run_query(filter_string=f"{column}:{filter_string}")

        for warning_substring in expected_warning:
            assert any(warning_substring.lower() in w.lower() for w in result.warnings)

        if expected_sample_ids is not None:
            sample_ids = result.get_unique_values("Sample_ID")
            assert sample_ids == expected_sample_ids or len(sample_ids) == len(expected_sample_ids)
        elif column == "Area" and filter_string == ">North":
            sample_ids = result.get_unique_values("Sample_ID")
            assert len(sample_ids) == 0


class TestSemicolonColumnTypeClassification:
    """Test that column type classification correctly handles semicolon-separated values.
    A column is numeric only if all parts of all semicolon-separated cells are valid numbers.
    If any part is non-numeric (e.g., '1-2' in '1;1-2'), the entire column is string."""

    def test_semicolon_cell_with_non_numeric_part_makes_column_string(
        self, sample_tsv_with_semicolon_mixed_type_column
    ):
        """Test that a column with a cell '1;1-2' should be treated as string because '1-2' is not a number."""
        manager = SidecarQueryManager(file=sample_tsv_with_semicolon_mixed_type_column)

        assert not pd.api.types.is_numeric_dtype(manager.df["Code"])
        assert not manager._is_semicolon_separated_numeric_column("Code")
        assert manager._is_mixed_type_column("Code")

    def test_semicolon_cell_with_non_numeric_part_warns_on_inequality(
        self, sample_tsv_with_semicolon_mixed_type_column
    ):
        """Test that inequality filter on a column broken by '1;1-2' should produce a warning."""
        manager = SidecarQueryManager(file=sample_tsv_with_semicolon_mixed_type_column)
        result = manager.run_query(filter_string="Code:>2")

        assert any("mixed types" in w.lower() and "Code" in w for w in result.warnings)
        assert any("comparison operators" in w.lower() for w in result.warnings)

    def test_semicolon_cell_with_non_numeric_part_string_matching_works(
        self, sample_tsv_with_semicolon_mixed_type_column
    ):
        """Test that string matching should still work on the mixed column. Filtering for '1-2' should matches cell value '1;1-2'."""
        manager = SidecarQueryManager(file=sample_tsv_with_semicolon_mixed_type_column)
        result = manager.run_query(filter_string="Code:1-2")

        sample_ids = result.get_unique_values("Sample_ID")
        assert "S1" in sample_ids
        assert len(sample_ids) == 1

    def test_semicolon_cell_with_non_numeric_part_single_numeric_match(
        self, sample_tsv_with_semicolon_mixed_type_column
    ):
        """Test that string matching for '3' on the mixed column should return S2 (exact string match)."""
        manager = SidecarQueryManager(file=sample_tsv_with_semicolon_mixed_type_column)
        result = manager.run_query(filter_string="Code:3")

        sample_ids = result.get_unique_values("Sample_ID")
        assert "S2" in sample_ids
        assert len(sample_ids) == 1

    def test_purely_numeric_semicolon_column_supports_numeric_ops(self, sample_tsv_with_semicolon_mixed_type_column):
        """Test that a column with only numeric semicolon values (e.g., '10;20;30') should support numeric operations."""
        manager = SidecarQueryManager(file=sample_tsv_with_semicolon_mixed_type_column)

        assert manager._is_semicolon_separated_numeric_column("PureNumericSemicolon")

        result = manager.run_query(filter_string="PureNumericSemicolon:>55")
        sample_ids = result.get_unique_values("Sample_ID")
        assert "S3" in sample_ids
        assert "S4" in sample_ids
        assert "S1" not in sample_ids
        assert "S2" not in sample_ids
        assert not any("string column" in w.lower() for w in result.warnings)

    def test_purely_numeric_semicolon_column_range_filter(self, sample_tsv_with_semicolon_mixed_type_column):
        """Test that a purely numeric semicolon column should support range operations."""
        manager = SidecarQueryManager(file=sample_tsv_with_semicolon_mixed_type_column)
        result = manager.run_query(filter_string="PureNumericSemicolon:25-45")

        sample_ids = result.get_unique_values("Sample_ID")
        assert "S1" in sample_ids
        assert "S2" in sample_ids
        assert "S3" not in sample_ids
        assert "S4" not in sample_ids
        assert not any("string column" in w.lower() for w in result.warnings)


class TestLoadFileValidation:
    """Test that SidecarQueryManager validates the same errors as the client-side
    MetadataTSVValidator in load_file()), before any queries are run.

    This ensures that even if a user skips the CLI validator, the server-side
    query engine catches the same formatting issues with clear error messages."""

    def test_commas_in_data_raises_at_tsv_load(self, tmp_path):
        """Test that commas in any cell value raises SidecarMetadataFormatError during load_file()."""
        tsv_content = "#Sample_ID\tArea\tWeight\nS1\tNorth,South\t12.5\nS2\tEast\t18.0\n"
        tsv_file = tmp_path / "commas.tsv"
        tsv_file.write_text(tsv_content)

        with pytest.raises(SidecarMetadataFormatError) as excinfo:
            SidecarQueryManager(file=tsv_file)
        assert "commas" in str(excinfo.value).lower()

    def test_commas_in_non_queried_column_raises_at_tsv_load(self, tmp_path):
        """Test that commas are caught in all columns when the tsv is loaded, not just the column being queried."""
        tsv_content = "#Sample_ID\tArea\tBadCol\nS1\tNorth\thas,comma\nS2\tEast\tclean\n"
        tsv_file = tmp_path / "commas_other_col.tsv"
        tsv_file.write_text(tsv_content)

        with pytest.raises(SidecarMetadataFormatError) as excinfo:
            SidecarQueryManager(file=tsv_file)
        assert "commas" in str(excinfo.value).lower()

    def test_duplicate_column_names_raises(self, tmp_path):
        """Test that duplicate column names raise SidecarMetadataFormatError during load_file().
        Without this check, pandas might silently rename them (e.g., 'Area', 'Area.1')."""
        tsv_content = "#Sample_ID\tArea\tArea\nS1\tNorth\tSouth\nS2\tEast\tWest\n"
        tsv_file = tmp_path / "duplicate_cols.tsv"
        tsv_file.write_text(tsv_content)

        with pytest.raises(SidecarMetadataFormatError) as excinfo:
            SidecarQueryManager(file=tsv_file)
        assert "duplicate" in str(excinfo.value).lower()

    def test_empty_column_name_raises(self, tmp_path):
        """Test that empty column names raise SidecarMetadataFormatError during load_file()."""
        tsv_content = "#Sample_ID\t\tWeight\nS1\tNorth\t12.5\nS2\tEast\t18.0\n"
        tsv_file = tmp_path / "empty_col.tsv"
        tsv_file.write_text(tsv_content)

        with pytest.raises(SidecarMetadataFormatError) as excinfo:
            SidecarQueryManager(file=tsv_file)
        assert "empty" in str(excinfo.value).lower()

    def test_semicolon_in_sample_id_raises(self, tmp_path):
        """Test that semicolons in Sample_ID values raise SidecarSampleIDError during load_file().
        Sample_ID must contain exactly one value per row."""
        tsv_content = "#Sample_ID\tArea\nS1;S2\tNorth\nS3\tEast\n"
        tsv_file = tmp_path / "semicolon_sample_id.tsv"
        tsv_file.write_text(tsv_content)

        with pytest.raises(SidecarSampleIDError) as excinfo:
            SidecarQueryManager(file=tsv_file)
        assert "semicolon" in str(excinfo.value).lower()

    def test_missing_sample_id_column_raises(self, sample_tsv_missing_sample_id_column):
        """Test that missing Sample_ID column raise SidecarColumnNotFoundError."""
        with pytest.raises(SidecarColumnNotFoundError):
            SidecarQueryManager(file=sample_tsv_missing_sample_id_column)

    def test_empty_sample_id_raises(self, sample_tsv_with_invalid_sample_ids):
        """Test that empty Sample_ID values raise SidecarSampleIDError."""
        with pytest.raises(SidecarSampleIDError):
            SidecarQueryManager(file=sample_tsv_with_invalid_sample_ids)

    def test_duplicate_sample_id_raises(self, sample_tsv_with_duplicate_sample_ids):
        """Test that duplicate Sample_ID values raise SidecarSampleIDError."""
        with pytest.raises(SidecarSampleIDError) as excinfo:
            SidecarQueryManager(file=sample_tsv_with_duplicate_sample_ids)
        assert "duplicate" in str(excinfo.value).lower()

    def test_valid_file_loads_successfully(self, sample_tsv_with_edge_cases):
        """
        Test that a TSV that follows DivBase requirements loads without errors.
        Use the edge case fixture to assert that these are fine as long as they all
        follow the DivBase requirements.
        """
        manager = SidecarQueryManager(file=sample_tsv_with_edge_cases)
        assert manager.df is not None
        assert "Sample_ID" in manager.df.columns
