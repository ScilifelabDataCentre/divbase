"""
Unit tests for SidecarQueryManager filtering
"""

import pytest

from divbase_api.services.queries import SidecarQueryManager


@pytest.fixture
def sample_tsv_with_numeric_data(tmp_path):
    """
    Create a temporary TSV file with numeric and string columns for testing.
    Includes semicolon-separated values in some cells.
    Note: Weight and Age columns have NO semicolons to ensure pandas infers them as numeric.
    Population column has semicolons but should still be numeric.
    """
    tsv_content = """#Sample_ID\tPopulation\tWeight\tAge\tArea
S1\t1\t20.0\t5.0\tNorth
S2\t2;4\t25.0\t10\tEast
S3\t3\t30.0\t15\tWest;South
S4\t4\t35.0\t20\tWest
S5\t5\t40.0\t25\tNorth
S6\t6\t45.0\t30\tEast
S7\t1;3;5\t50.0\t35\tSouth
S8\t2\t55.0\t40\tWest
S9\t7\t62.0\t45\tNorth
S10\t8\t70.0\t52\tEast
"""
    tsv_file = tmp_path / "test_metadata.tsv"
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
        assert "S8" in sample_ids  # Weight: 55.0
        assert "S9" in sample_ids  # Weight: 62.0
        assert "S10" in sample_ids  # Weight: 70.0

    def test_greater_than_or_equal(self, sample_tsv_with_numeric_data):
        """Test >= operator returns correct samples."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Weight:>=50")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 4
        assert "S7" in sample_ids  # Weight: 50.0
        assert "S8" in sample_ids  # Weight: 55.0
        assert "S9" in sample_ids  # Weight: 62.0
        assert "S10" in sample_ids  # Weight: 70.0

    def test_less_than(self, sample_tsv_with_numeric_data):
        """Test < operator returns correct samples."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Age:<15")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 2
        assert "S1" in sample_ids  # Age: 5
        assert "S2" in sample_ids  # Age: 10

    def test_less_than_or_equal(self, sample_tsv_with_numeric_data):
        """Test <= operator returns correct samples."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Age:<=15")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 3
        assert "S1" in sample_ids  # Age: 5
        assert "S2" in sample_ids  # Age: 10
        assert "S3" in sample_ids  # Age: 15

    def test_inequality_on_weight_column(self, sample_tsv_with_numeric_data):
        """Test inequality on Weight column (no semicolons, pure numeric)."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Weight:>60")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 2
        assert "S9" in sample_ids  # Weight: 62.0
        assert "S10" in sample_ids  # Weight: 70.0

    def test_inequality_on_age_column(self, sample_tsv_with_numeric_data):
        """Test inequality on Age column (no semicolons, pure numeric)."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Age:>=40")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 3
        assert "S8" in sample_ids  # Age: 40.0
        assert "S9" in sample_ids  # Age: 45.0
        assert "S10" in sample_ids  # Age: 52.0


class TestNumericalFilteringRanges:
    """Test range filtering on numerical columns."""

    def test_simple_range(self, sample_tsv_with_numeric_data):
        """Test inclusive range filtering."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Weight:30-45")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 4
        assert "S3" in sample_ids  # Weight: 30
        assert "S4" in sample_ids  # Weight: 35
        assert "S5" in sample_ids  # Weight: 40
        assert "S6" in sample_ids  # Weight: 45

    def test_range_boundaries_inclusive(self, sample_tsv_with_numeric_data):
        """Test that range boundaries are inclusive."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Age:20-30")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 3
        assert "S4" in sample_ids  # Age: 20 (lower boundary)
        assert "S5" in sample_ids  # Age: 25
        assert "S6" in sample_ids  # Age: 30 (upper boundary)

    def test_range_on_weight_column(self, sample_tsv_with_numeric_data):
        """Test range filtering on Weight column (no semicolons, pure numeric)."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Weight:40-60")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 4
        assert "S5" in sample_ids  # Weight: 40.0
        assert "S6" in sample_ids  # Weight: 45.0
        assert "S7" in sample_ids  # Weight: 50.0
        assert "S8" in sample_ids  # Weight: 55.0

    def test_narrow_range(self, sample_tsv_with_numeric_data):
        """Test a narrow range returns only matching samples."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Age:20-30")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 3
        assert "S4" in sample_ids  # Age: 20.0
        assert "S5" in sample_ids  # Age: 25.0
        assert "S6" in sample_ids  # Age: 30.0


class TestNumericalFilteringDiscreteValues:
    """Test discrete value filtering on numerical columns."""

    def test_single_discrete_value(self, sample_tsv_with_numeric_data):
        """Test filtering with a single discrete value."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Weight:50")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 1
        assert "S7" in sample_ids  # Weight: 50

    def test_multiple_discrete_values(self, sample_tsv_with_numeric_data):
        """Test filtering with multiple discrete values (OR logic)."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Weight:20,30,50")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 3
        assert "S1" in sample_ids  # Weight: 20
        assert "S3" in sample_ids  # Weight: 30
        assert "S7" in sample_ids  # Weight: 50

    def test_discrete_values_with_semicolon_separated_cells(self, sample_tsv_with_numeric_data):
        """Test discrete value filtering on Population column (string column with semicolons)."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Population:1,3,5")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 4
        assert "S1" in sample_ids  # Population: 1
        assert "S3" in sample_ids  # Population: 3
        assert "S5" in sample_ids  # Population: 5
        assert "S7" in sample_ids  # Population: 1;3;5 (string matches "1", "3", and "5")

    def test_discrete_values_match_any_semicolon_value(self, sample_tsv_with_numeric_data):
        """Test that discrete filtering matches if ANY semicolon value matches (string matching on Population)."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Population:4")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 2
        assert "S2" in sample_ids  # Population: 2;4 (string "4" matches)
        assert "S4" in sample_ids  # Population: 4

    def test_discrete_age_values(self, sample_tsv_with_numeric_data):
        """Test discrete filtering on Age column."""
        manager = SidecarQueryManager(file=sample_tsv_with_numeric_data)
        result = manager.run_query(filter_string="Age:10,25,40")

        sample_ids = result.get_unique_values("Sample_ID")
        assert len(sample_ids) == 3
        assert "S2" in sample_ids  # Age: 10
        assert "S5" in sample_ids  # Age: 25
        assert "S8" in sample_ids  # Age: 40
