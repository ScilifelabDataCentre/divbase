"""
Unit tests for SidecarQueryManager filtering.

Fixtures use Python list notation for multi-value cells (e.g. [2, 4] instead
of the old semicolon-separated format 2;4).  Shared fixtures live in
tests/unit/conftest.py.
"""

import pandas as pd
import pytest

from divbase_api.services import queries as queries_module
from divbase_api.services.queries import SidecarQueryManager
from divbase_lib.exceptions import (
    SidecarColumnNotFoundError,
    SidecarMetadataFormatError,
    SidecarNoDataLoadedError,
    SidecarSampleIDError,
)
from divbase_lib.metadata_validator import (
    MetadataValidationResult,
    ValidationCategory,
    ValidationMessage,
)


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
        """Test inequality on Weight column (pure numeric scalars)."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Weight:>60")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 2
        assert "S9" in sample_ids
        assert "S10" in sample_ids

    def test_inequality_on_age_column(self, sample_tsv_with_numeric_data):
        """Test inequality on Age column (pure numeric scalars)."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Age:>=40")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 3
        assert "S8" in sample_ids
        assert "S9" in sample_ids
        assert "S10" in sample_ids

    def test_inequality_on_single_numeric_value_column(self, sample_tsv_with_numeric_data):
        """Test inequality on single-value numeric column (no list cells)."""
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
        """Test range filtering on Weight column (pure numeric scalars)."""
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

    def test_discrete_values_on_list_column(self, sample_tsv_with_numeric_data):
        """Test discrete value filtering on Population column (has list cells like [2, 4])."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Population:1,3,5")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 4
        assert "S1" in sample_ids
        assert "S3" in sample_ids
        assert "S5" in sample_ids
        assert "S7" in sample_ids

    def test_discrete_values_match_any_list_element(self, sample_tsv_with_numeric_data):
        """Test that discrete filtering matches if ANY list element matches."""
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


class TestListValueNumericFiltering:
    """Test that inequalities and ranges work on columns with Python list multi-value cells."""

    def test_inequality_on_list_column(self, sample_tsv_with_numeric_data):
        """Test that > operator works on Population column (has list cells like [2, 4])."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Population:>4")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 5
        assert "S5" in sample_ids
        assert "S6" in sample_ids
        assert "S7" in sample_ids
        assert "S9" in sample_ids
        assert "S10" in sample_ids

    def test_inequality_less_than_on_list_column(self, sample_tsv_with_numeric_data):
        """Test that < operator works on Population column."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Population:<3")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 4
        assert "S1" in sample_ids
        assert "S2" in sample_ids
        assert "S7" in sample_ids
        assert "S8" in sample_ids

    def test_range_on_list_column(self, sample_tsv_with_numeric_data):
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

    def test_combined_inequality_and_discrete_on_list_column(self, sample_tsv_with_numeric_data):
        """Test combining inequality and discrete values on Population column."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Population:>6,2")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 4
        assert "S2" in sample_ids
        assert "S8" in sample_ids
        assert "S9" in sample_ids
        assert "S10" in sample_ids

    def test_range_with_list_values_at_boundaries(self, sample_tsv_with_numeric_data):
        """Test that range boundaries work correctly with list cells."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Population:1-3")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 5
        assert "S1" in sample_ids
        assert "S2" in sample_ids
        assert "S3" in sample_ids
        assert "S7" in sample_ids
        assert "S8" in sample_ids

    def test_not_operator_with_list_values(self, sample_tsv_with_numeric_data):
        """Test NOT operator (!) with list-value column: Population:>3,!5 should exclude 5."""
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
    """Test string column filtering with single and list multi-value cells."""

    def test_single_string_value_column(self, sample_tsv_with_numeric_data):
        """Test filtering on a string column with single values."""
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

    def test_list_string_column(self, sample_tsv_with_numeric_data):
        """Test filtering on string column with list cells (e.g. ["West", "South", "East"])."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Area:West")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 3
        assert "S4" in sample_ids
        assert "S8" in sample_ids
        assert "S3" in sample_ids

    def test_list_string_column_any_match(self, sample_tsv_with_numeric_data):
        """Test that filtering matches if ANY list element matches."""
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
        """Test that a column with mixed numeric and non-numeric scalar values is treated as a string column.
        The MixedTypes column has values like 1, 'two', 3, 'string4'.
        When treated as string, filtering for '1' should match cells containing '1'."""
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
        """Test that filtering for Unicode values like 'Göteborg' works and returns correct samples."""
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

    def test_list_string_filtering_matches_element(self, sample_tsv_with_edge_cases):
        """Test that filtering on a string column with list cells matches individual list elements."""
        manager = SidecarQueryManager(file=sample_tsv_with_edge_cases)
        result = manager.run_query(filter_string="PureStrings:South")
        sample_ids = result.get_unique_values("Sample_ID")
        assert "S1" in sample_ids
        assert "S3" in sample_ids


class TestSampleIDValidation:
    """Test Sample_ID validation during file loading."""

    def test_empty_sample_id_raises_error(self, sample_tsv_with_invalid_sample_ids):
        """Test that empty Sample_ID values raise SidecarSampleIDError directly during file load."""
        with pytest.raises(SidecarSampleIDError) as excinfo:
            SidecarQueryManager(file=sample_tsv_with_invalid_sample_ids)
        assert "sample_id" in str(excinfo.value).lower() and "empty" in str(excinfo.value).lower()

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
        assert "Sample_ID" in str(excinfo.value)


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

    def test_negative_numbers_in_list_cells(self, sample_tsv_with_numeric_data):
        """Test that negative numbers in list cells (e.g. [-3.5, -2.1]) work correctly."""
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


class TestListColumnTypeClassification:
    """Test that column type classification correctly handles list and mixed-type cells.
    A column is numeric only if all cells (scalars and list elements) are valid numbers.
    If any cell is non-numeric, the column is treated as string/mixed-type."""

    def test_non_numeric_scalar_makes_column_mixed(self, sample_tsv_with_list_mixed_type_column):
        """Test that a column with a non-numeric scalar ('1-2') among numeric values is mixed-type."""
        manager = SidecarQueryManager(file=sample_tsv_with_list_mixed_type_column)

        assert "Code" in manager.mixed_type_columns

    def test_mixed_column_warns_on_inequality(self, sample_tsv_with_list_mixed_type_column):
        """Test that inequality filter on a mixed-type column should produce a warning."""
        manager = SidecarQueryManager(file=sample_tsv_with_list_mixed_type_column)
        result = manager.run_query(filter_string="Code:>2")

        assert any("mixed types" in w.lower() and "Code" in w for w in result.warnings)
        assert any("comparison operators" in w.lower() for w in result.warnings)

    def test_mixed_column_string_matching_works(self, sample_tsv_with_list_mixed_type_column):
        """Test that string matching works on a mixed-type column. Filtering for '1-2' should match S1."""
        manager = SidecarQueryManager(file=sample_tsv_with_list_mixed_type_column)
        result = manager.run_query(filter_string="Code:1-2")

        sample_ids = result.get_unique_values("Sample_ID")
        assert "S1" in sample_ids
        assert len(sample_ids) == 1

    def test_mixed_column_single_numeric_match(self, sample_tsv_with_list_mixed_type_column):
        """Test that string matching for '3' on the mixed column should return S2 (exact string match)."""
        manager = SidecarQueryManager(file=sample_tsv_with_list_mixed_type_column)
        result = manager.run_query(filter_string="Code:3")

        sample_ids = result.get_unique_values("Sample_ID")
        assert "S2" in sample_ids
        assert len(sample_ids) == 1

    def test_purely_numeric_list_column_supports_numeric_ops(self, sample_tsv_with_list_mixed_type_column):
        """Test that a column with numeric list cells (e.g. [10, 20, 30]) supports numeric operations."""
        manager = SidecarQueryManager(file=sample_tsv_with_list_mixed_type_column)

        assert "PureNumericList" in manager.numeric_columns

        result = manager.run_query(filter_string="PureNumericList:>55")
        sample_ids = result.get_unique_values("Sample_ID")
        assert "S3" in sample_ids
        assert "S4" in sample_ids
        assert "S1" not in sample_ids
        assert "S2" not in sample_ids
        assert not any("string column" in w.lower() for w in result.warnings)

    def test_purely_numeric_list_column_range_filter(self, sample_tsv_with_list_mixed_type_column):
        """Test that a purely numeric list column supports range operations."""
        manager = SidecarQueryManager(file=sample_tsv_with_list_mixed_type_column)
        result = manager.run_query(filter_string="PureNumericList:25-45")

        sample_ids = result.get_unique_values("Sample_ID")
        assert "S1" in sample_ids
        assert "S2" in sample_ids
        assert "S3" not in sample_ids
        assert "S4" not in sample_ids
        assert not any("string column" in w.lower() for w in result.warnings)


class TestLoadFileValidation:
    """Test that SidecarQueryManager validates the same errors as the client-side
    ClientSideClientSideMetadataTSVValidator in load_file()), before any queries are run.

    This ensures that even if a user skips the CLI validator, the server-side
    query engine catches the same formatting issues with clear error messages."""

    def test_commas_in_mixed_numeric_column_detected_during_query(self, tmp_path):
        """Test that commas in a column with mixed numeric and non-numeric values trigger mixed-type warning."""
        tsv_content = "#Sample_ID\tPopulation\tWeight\nS1\t1,2\t12.5\nS2\t5\t18.0\n"
        tsv_file = tmp_path / "commas.tsv"
        tsv_file.write_text(tsv_content)

        manager = SidecarQueryManager(file=tsv_file)
        result = manager.run_query("Population:1,2")
        assert any("mixed types" in w.lower() for w in result.warnings)

    def test_commas_in_pure_string_column_no_mixed_type_warning(self, tmp_path):
        """Test that commas in a pure string column don't trigger mixed-type warnings."""
        tsv_content = "#Sample_ID\tCode\nS1\t1,2\nS2\t3,4\nS3\t5,6\n"
        tsv_file = tmp_path / "commas.tsv"
        tsv_file.write_text(tsv_content)

        manager = SidecarQueryManager(file=tsv_file)
        result = manager.run_query("Code:1,2")

        assert not any("mixed types" in w.lower() for w in result.warnings)

    def test_commas_in_cells_produce_format_warning(self, tmp_path):
        """Test that commas in plain string cells produce FORMAT warnings during load."""
        tsv_content = "#Sample_ID\tCode\nS1\t1,2\nS2\t3,4\nS3\t5,6\n"
        tsv_file = tmp_path / "commas.tsv"
        tsv_file.write_text(tsv_content)

        manager = SidecarQueryManager(file=tsv_file)
        assert any("comma" in w.lower() for w in manager.warnings)

    def test_duplicate_column_names_raises(self, tmp_path):
        """Test that duplicate column names raise SidecarMetadataFormatError during load_file().
        Without this check, pandas might silently rename them (e.g., 'Area', 'Area.1')."""
        tsv_content = "#Sample_ID\tArea\tArea\nS1\tNorth\tSouth\nS2\tEast\tWest\n"
        tsv_file = tmp_path / "duplicate_cols.tsv"
        tsv_file.write_text(tsv_content)

        with pytest.raises(SidecarMetadataFormatError) as excinfo:
            SidecarQueryManager(file=tsv_file)
        assert "duplicate" in str(excinfo.value).lower()

    def test_duplicate_column_names_after_stripping_hash_raises(self, tmp_path):
        """Test that duplicate column names are caught even when one has '#' and one doesn't.
        For example, '#Sample_ID' and 'Sample_ID' should be detected as duplicates."""
        tsv_content = "#Sample_ID\tSample_ID\tPopulation\nS1\tS1_dup\t1\nS2\tS2_dup\t2\n"
        tsv_file = tmp_path / "duplicate_sample_id_cols.tsv"
        tsv_file.write_text(tsv_content)

        with pytest.raises(SidecarMetadataFormatError) as excinfo:
            SidecarQueryManager(file=tsv_file)
        assert "duplicate" in str(excinfo.value).lower()
        assert "Sample_ID" in str(excinfo.value)

    def test_empty_column_name_raises(self, tmp_path):
        """Test that empty column names raise SidecarMetadataFormatError during load_file()."""
        tsv_content = "#Sample_ID\t\tWeight\nS1\tNorth\t12.5\nS2\tEast\t18.0\n"
        tsv_file = tmp_path / "empty_col.tsv"
        tsv_file.write_text(tsv_content)

        with pytest.raises(SidecarMetadataFormatError) as excinfo:
            SidecarQueryManager(file=tsv_file)
        assert "empty" in str(excinfo.value).lower()

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

    def test_invalid_list_syntax_raises(self, tmp_path):
        """Test that invalid Python list syntax in cells raises an error during load."""
        tsv_content = "#Sample_ID\tScores\nS1\t[1, 2\nS2\t5\n"
        tsv_file = tmp_path / "bad_list.tsv"
        tsv_file.write_text(tsv_content)

        with pytest.raises(SidecarMetadataFormatError) as excinfo:
            SidecarQueryManager(file=tsv_file)
        assert "invalid" in str(excinfo.value).lower()


class TestQuotedFilterValues:
    """Test that quoted filter values allow querying for strings containing commas."""

    def test_quoted_filter_value_with_comma(self, tmp_path):
        """Test that a quoted filter value like '"North,South"' matches the literal string."""
        tsv_content = '#Sample_ID\tArea\nS1\t"North,South"\nS2\tEast\nS3\tNorth\n'
        tsv_file = tmp_path / "quoted.tsv"
        tsv_file.write_text(tsv_content)

        manager = SidecarQueryManager(file=tsv_file)
        result = manager.run_query(filter_string='Area:"North,South"')
        sample_ids = result.get_unique_values("Sample_ID")
        assert "S1" in sample_ids
        assert len(sample_ids) == 1

    def test_negation_with_quoted_value(self, tmp_path):
        """Test that negation works with quoted values: Area:!"North,South" excludes the literal."""
        tsv_content = '#Sample_ID\tArea\nS1\t"North,South"\nS2\tEast\nS3\tNorth\n'
        tsv_file = tmp_path / "quoted_neg.tsv"
        tsv_file.write_text(tsv_content)

        manager = SidecarQueryManager(file=tsv_file)
        result = manager.run_query(filter_string='Area:!"North,South"')
        sample_ids = result.get_unique_values("Sample_ID")
        assert "S2" in sample_ids
        assert "S3" in sample_ids
        assert "S1" not in sample_ids

    def test_unquoted_comma_splits_filter_values(self, tmp_path):
        """Test that unquoted commas still split filter values as OR conditions."""
        tsv_content = "#Sample_ID\tArea\nS1\tNorth\nS2\tEast\nS3\tSouth\n"
        tsv_file = tmp_path / "unquoted.tsv"
        tsv_file.write_text(tsv_content)

        manager = SidecarQueryManager(file=tsv_file)
        result = manager.run_query(filter_string="Area:North,East")
        sample_ids = result.get_unique_values("Sample_ID")
        assert "S1" in sample_ids
        assert "S2" in sample_ids
        assert "S3" not in sample_ids

    def test_quoted_filter_value_with_comma_and_space(self, tmp_path):
        """Test that a quoted filter value can match a literal value like 'North, South'."""
        tsv_content = '#Sample_ID\tArea\nS1\t"North, South"\nS2\tEast\nS3\tNorth\n'
        tsv_file = tmp_path / "quoted_comma_space.tsv"
        tsv_file.write_text(tsv_content)

        manager = SidecarQueryManager(file=tsv_file)
        result = manager.run_query(filter_string='Area:"North, South"')
        sample_ids = result.get_unique_values("Sample_ID")
        assert sample_ids == ["S1"]

    def test_unquoted_filter_with_comma_and_space_does_not_match_literal(self, tmp_path):
        """Test that unquoted comma-separated filter values do not match a literal 'North, South' cell."""
        tsv_content = '#Sample_ID\tArea\nS1\t"North, South"\nS2\tNorth\nS3\tSouth\n'
        tsv_file = tmp_path / "unquoted_comma_space.tsv"
        tsv_file.write_text(tsv_content)

        manager = SidecarQueryManager(file=tsv_file)
        result = manager.run_query(filter_string="Area:North, South")
        sample_ids = sorted(result.get_unique_values("Sample_ID"))
        assert sample_ids == ["S2", "S3"]
        assert "S1" not in sample_ids

    def test_embedded_double_quote_value_is_not_queryable_with_current_parser(self, tmp_path):
        """Document current behavior: values containing literal double quotes are not queryable."""
        tsv_content = '#Sample_ID\tArea\nS1\t"He said ""Hi"""\nS2\tEast\n'
        tsv_file = tmp_path / "embedded_quotes.tsv"
        tsv_file.write_text(tsv_content)

        manager = SidecarQueryManager(file=tsv_file)
        result = manager.run_query(filter_string='Area:"He said ""Hi"""')
        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 0
        assert any("No results for the filter" in w and "column 'Area'" in w for w in result.warnings)


class TestValidationCategoryDispatch:
    """Test how SidecarQueryManager maps SharedMetadataValidator categories to behavior."""

    @staticmethod
    def _patch_shared_validator(monkeypatch, result):
        class FakeSharedMetadataValidator:
            def __init__(self, file_path, project_samples=None, skip_dimensions_check=False):
                self.file_path = file_path
                self.project_samples = project_samples
                self.skip_dimensions_check = skip_dimensions_check

            def load_and_validate(self):
                return result

        monkeypatch.setattr(queries_module, "SharedMetadataValidator", FakeSharedMetadataValidator)

    def test_file_read_category_raises_no_data_loaded(self, monkeypatch, tmp_path):
        result = MetadataValidationResult(
            errors=[ValidationMessage(ValidationCategory.FILE_READ, "failed to read")],
        )
        self._patch_shared_validator(monkeypatch, result)
        with pytest.raises(SidecarNoDataLoadedError):
            SidecarQueryManager(file=tmp_path / "dummy.tsv")

    def test_sample_id_column_category_raises_column_not_found(self, monkeypatch, tmp_path):
        result = MetadataValidationResult(
            errors=[ValidationMessage(ValidationCategory.SAMPLE_ID_COLUMN, "Sample_ID column is required")],
        )
        self._patch_shared_validator(monkeypatch, result)
        with pytest.raises(SidecarColumnNotFoundError):
            SidecarQueryManager(file=tmp_path / "dummy.tsv")

    def test_sample_id_value_category_raises_sample_id_error(self, monkeypatch, tmp_path):
        result = MetadataValidationResult(
            errors=[ValidationMessage(ValidationCategory.SAMPLE_ID_VALUE, "Sample_ID is empty")],
        )
        self._patch_shared_validator(monkeypatch, result)
        with pytest.raises(SidecarSampleIDError):
            SidecarQueryManager(file=tmp_path / "dummy.tsv")

    def test_other_error_category_raises_metadata_format_error(self, monkeypatch, tmp_path):
        result = MetadataValidationResult(
            errors=[ValidationMessage(ValidationCategory.HEADER, "Duplicate column names found")],
        )
        self._patch_shared_validator(monkeypatch, result)
        with pytest.raises(SidecarMetadataFormatError):
            SidecarQueryManager(file=tmp_path / "dummy.tsv")

    def test_only_dimensions_and_format_warnings_are_forwarded(self, monkeypatch, tmp_path):
        result = MetadataValidationResult(
            warnings=[
                ValidationMessage(ValidationCategory.DIMENSIONS, "dimension warning"),
                ValidationMessage(ValidationCategory.FORMAT, "format warning"),
                ValidationMessage(ValidationCategory.TYPE_CLASSIFICATION, "type warning"),
            ],
            df=pd.DataFrame({"Sample_ID": ["S1"]}),
        )
        self._patch_shared_validator(monkeypatch, result)
        manager = SidecarQueryManager(file=tmp_path / "dummy.tsv")
        assert "dimension warning" in manager.warnings
        assert "format warning" in manager.warnings
        assert "type warning" not in manager.warnings
