"""
Unit tests for SidecarQueryManager filtering
"""

import pytest

from divbase_api.services.queries import SidecarQueryManager
from divbase_lib.exceptions import SidecarColumnNotFoundError, SidecarInvalidFilterError, SidecarSampleIDError


@pytest.fixture
def sample_tsv_with_numeric_data(tmp_path):
    """
    Create a temporary TSV file with numeric and string columns for testing.
    Includes semicolon-separated values in some cells. Includes both int and float
    numeric values to test that both are detected as numeric.
    """
    # Keep indentation like this to ensure that leading spaces in column 1 does not cause issues.
    tsv_content = """#Sample_ID\tPopulation\tWeight\tAge\tArea\tSingleNumber\tSingleString
S1\t1\t20.2\t5.0\tNorth\t100\tString
S2\t2;4\t25.0\t10\tEast\t200\tStrings
S3\t3\t30.8\t15\tWest;South;East\t300\tSting
S4\t4\t35.1\t20\tWest\t400\tStings
S5\t5\t40.0\t25\tNorth\t500\tThing
S6\t6\t45.4\t30\tEast\t600\tThings
S7\t1;3;5\t50.9\t35\tSouth\t700\tStrong
S8\t2\t55.2\t40\tWest\t800\tStrung
S9\t7\t62.6\t45\tNorth\t900\tStang
S10\t8\t70.7\t52\tEast\t1000\tSong
"""
    tsv_file = tmp_path / "test_metadata.tsv"
    tsv_file.write_text(tsv_content)
    return tsv_file


@pytest.fixture
def sample_tsv_with_edge_cases(tmp_path):
    """
    Create a temporary TSV file to test edge cases:
    1. "string;string;string" - OK (pure strings)
    2. "1;two;5" - FAIL (mixed numeric and non-numeric should raise exception)
    3. String values containing numbers like "1string" - OK (inferred as string)
    4. Column with commas should raise SidecarInvalidFilterError

    Commas are NOT allowed in divbase TSV format.
    Note that S2 and S3 have leading/trailing whitespace in the Sample_ID and the code should handle that by stripping whitespace.
    """
    tsv_content = """#Sample_ID\tPureStrings\tMixedTypes\tSingleString\tSingleNumber\tUnicodeStrings\tWithCommas\tStringWithHyphen\tNumericalWithHyphen
S1\tNorth;South;East\t1;two;5\tWest\t100\tStockholm;Göteborg\tNorth,South\tNorth-East\t1-2
S2 \tWest;East;North\t2;three;6\tNorth\t200\tMalmö;Uppsala\tWest,East\tSouth-West\t2-3
 S3\tSouth\t3\tEast\t300\tKöpenhamn;København\tNorth,\tNorth-North-West\t3-4
S4\t1string\tstring4\tString5\t400\tHumlebæk\t,South\tEast-South-East\t4-5
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

    def test_column_with_commas_raises(self, sample_tsv_with_edge_cases):
        """Test that a column containing commas raises SidecarInvalidFilterError."""
        manager = SidecarQueryManager(file=sample_tsv_with_edge_cases)
        with pytest.raises(SidecarInvalidFilterError) as excinfo:
            manager.run_query(filter_string="WithCommas:foo")
        assert "contains commas" in str(excinfo.value)

    def test_mixed_types_should_fail(self, sample_tsv_with_edge_cases):
        """Test that a column with mixed numeric and non-numeric values raises an error."""
        manager = SidecarQueryManager(file=sample_tsv_with_edge_cases)
        with pytest.raises(SidecarInvalidFilterError):
            manager.run_query(filter_string="MixedTypes:1")

    def test_strings_with_commas(self, sample_tsv_with_edge_cases):
        """Test that a column with strings containing commas is correctly handled."""
        manager = SidecarQueryManager(file=sample_tsv_with_edge_cases)
        result = manager.run_query(filter_string="StringsWithCommas:Region1")
        sample_ids = result.get_unique_values("Sample_ID")
        assert "S1" in sample_ids

    def test_numbers_with_comma(self, sample_tsv_with_edge_cases):
        """Test that a column with numeric values containing commas is treated as string type (since comma is not a numeric character)."""
        manager = SidecarQueryManager(file=sample_tsv_with_edge_cases)
        result = manager.run_query(filter_string="NumbersWithComma:1")
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

    def test_multi_column_single_string_and_mixed_types_should_fail(self, sample_tsv_with_edge_cases):
        """Test that filtering on SingleString (valid) and MixedTypes (invalid) will fail due to mixed types."""
        manager = SidecarQueryManager(file=sample_tsv_with_edge_cases)
        with pytest.raises(SidecarInvalidFilterError):
            manager.run_query(filter_string="SingleString:String5;MixedTypes:string4")

    def test_multi_column_single_number_and_mixed_types_should_fail(self, sample_tsv_with_edge_cases):
        """Test that filtering on SingleNumber (valid) and MixedTypes (invalid) will fail due to mixed types."""
        manager = SidecarQueryManager(file=sample_tsv_with_edge_cases)
        with pytest.raises(SidecarInvalidFilterError):
            manager.run_query(filter_string="SingleNumber:400;MixedTypes:string4")

    def test_multi_column_single_string_and_pure_strings(self, sample_tsv_with_edge_cases):
        """Test that filtering on SingleString and PureStrings (both valid string columns) will pass."""
        manager = SidecarQueryManager(file=sample_tsv_with_edge_cases)
        result = manager.run_query(filter_string="SingleString:String5;PureStrings:1string")
        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 1
        assert "S4" in sample_ids

    def test_multi_column_single_number_and_numbers_with_comma(self, sample_tsv_with_edge_cases):
        """Test that filtering on SingleNumber (numeric) and NumbersWithComma (treated as string due to comma) will pass."""
        manager = SidecarQueryManager(file=sample_tsv_with_edge_cases)
        result = manager.run_query(filter_string="SingleNumber:400;NumbersWithComma:string3")
        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 1
        assert "S4" in sample_ids

    def test_multi_column_with_unicode(self, sample_tsv_with_edge_cases):
        """Test that multi-column filtering works with unicode strings, but raises error if commas are present."""
        manager = SidecarQueryManager(file=sample_tsv_with_edge_cases)
        with pytest.raises(SidecarInvalidFilterError):
            manager.run_query(filter_string="UnicodeStrings:København;WithCommas:North")

    def test_hyphens_allowed_in_string_values(self, sample_tsv_with_edge_cases):
        """Test that hyphens are allowed in string values and can be filtered correctly."""
        manager = SidecarQueryManager(file=sample_tsv_with_edge_cases)
        result = manager.run_query(filter_string="StringWithHyphen:South-West")
        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 1
        assert "S2" in sample_ids

    def test_hyphens_in_numerical_column_raises(self, sample_tsv_with_edge_cases):
        """Test that hyphens are allowed in string columns but not in numerical columns."""
        manager = SidecarQueryManager(file=sample_tsv_with_edge_cases)
        with pytest.raises(SidecarInvalidFilterError) as excinfo:
            manager.run_query(filter_string="NumericalWithHyphen:2")
        assert "Column 'NumericalWithHyphen' contains value '1-2' with a hyphen at row 0." in str(excinfo.value)

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
