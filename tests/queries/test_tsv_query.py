import logging

import pytest
from conftest import create_sidecar_manager

from divbase_tools.exceptions import (
    SidecarInvalidFilterError,
)
from divbase_tools.queries import SidecarQueryManager


@pytest.mark.unit
def test_sidecar_query_manager_init(sample_tsv_file):
    """Test that SidecarQueryManager initializes properly and loads file data."""
    manager = SidecarQueryManager(file=sample_tsv_file)

    assert manager.df is not None, "DataFrame should be populated upon initialization"
    assert len(manager.df) > 0, "DataFrame should contain data"
    assert "Sample_ID" in manager.df.columns, "Expected columns should be present"


@pytest.mark.unit
def test_sidecar_query_manager_chain(sample_tsv_file):
    """Test that method chaining works properly."""
    query_result = SidecarQueryManager(file=sample_tsv_file).run_query(filter_string="Population:Pop1").query_result

    assert len(query_result) == 2, "Method chaining should work correctly"


@pytest.mark.unit
def test_tsv_query_empty_filter(sample_tsv_file, caplog):
    """Test that empty filter returns all records with a warning."""
    with caplog.at_level(logging.WARNING):
        manager = create_sidecar_manager(sample_tsv_file)
        manager.run_query(filter_string="")
        query_result = manager.query_result
        query_message = manager.query_message

    assert len(query_result) == 5, "Empty filter should return all records"
    assert "Empty filter provided - returning ALL records" in caplog.text
    assert "ALL RECORDS" in query_message


@pytest.mark.unit
def test_tsv_query_none_filter(sample_tsv_file):
    """Test that None filter raises an appropriate error."""
    manager = create_sidecar_manager(sample_tsv_file)
    with pytest.raises(SidecarInvalidFilterError) as exc_info:
        manager.run_query(filter_string=None)

    assert "Filter cannot be None" in str(exc_info.value)


@pytest.mark.unit
def test_tsv_query_column_not_found(sample_tsv_file, caplog):
    """Test handling of filter with non-existent column."""
    with caplog.at_level(logging.WARNING):
        manager = create_sidecar_manager(sample_tsv_file)
        manager.run_query(filter_string="NonExistentColumn:Value")
        query_result = manager.query_result
        query_message = manager.query_message

    assert len(query_result) == 5, "Should return all records when column not found"
    assert "Column 'NonExistentColumn' not found in the TSV file" in caplog.text
    assert query_message == "Invalid filter conditions - returning ALL records"


@pytest.mark.unit
def test_tsv_query_value_not_found(sample_tsv_file, caplog):
    """Test handling of filter with non-existent values."""
    with caplog.at_level(logging.WARNING):
        manager = create_sidecar_manager(sample_tsv_file)
        manager.run_query(filter_string="Population:NonExistentPop")
        query_result = manager.query_result

    assert len(query_result) == 0, "Should return empty DataFrame when no values match"
    assert "None of the values ['NonExistentPop'] were found in column 'Population'" in caplog.text


@pytest.mark.unit
def test_tsv_query_invalid_filter_format(sample_tsv_file, caplog):
    """Test handling of incorrectly formatted filter."""
    with caplog.at_level(logging.WARNING):
        manager = create_sidecar_manager(sample_tsv_file)
        with pytest.raises(SidecarInvalidFilterError) as exc_info:
            manager.run_query(filter_string="InvalidFormatNoColon")
        assert "Invalid filter format" in str(exc_info.value)


@pytest.mark.unit
def test_tsv_query_filename_filter(sample_tsv_file):
    """Test filtering by Filename column."""
    manager = create_sidecar_manager(sample_tsv_file)
    manager.run_query(filter_string="Filename:file1.vcf.gz")
    query_result = manager.query_result

    assert len(query_result) == 2, "Should return records with matching filename"
    assert list(query_result["Sample_ID"]) == ["S1", "S2"], "Should return correct samples"


@pytest.mark.unit
def test_tsv_query_multiple_conditions(sample_tsv_file):
    """Test filtering with multiple conditions (intersect query)."""
    manager = create_sidecar_manager(sample_tsv_file)
    manager.run_query(filter_string="Population:Pop1;Sex:F")
    query_result = manager.query_result

    assert len(query_result) == 1, "Should return records matching both conditions"
    assert list(query_result["Sample_ID"]) == ["S2"], "Should return correct sample"


@pytest.mark.unit
def test_tsv_query_multiple_values(sample_tsv_file):
    """Test filtering with multiple values for a single column."""
    manager = create_sidecar_manager(sample_tsv_file)
    manager.run_query(filter_string="Population:Pop1,Pop3")
    query_result = manager.query_result

    assert len(query_result) == 3, "Should return records matching any value"
    assert sorted(list(query_result["Sample_ID"])) == ["S1", "S2", "S5"], "Should return correct samples"


@pytest.mark.unit
def test_tsv_query_strip_hash_headers(tsv_with_hash_headers):
    """Test that column headers with # prefix are handled correctly."""
    manager = create_sidecar_manager(tsv_with_hash_headers)
    manager.run_query(filter_string="Sample_ID:S1")
    query_result = manager.query_result

    assert len(query_result) == 1, "Should strip # from headers and match correctly"
    assert list(query_result.columns) == ["Sample_ID", "Population", "Filename"], "Should strip # from headers"
