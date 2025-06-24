import logging

import pytest

from divbase_tools.queries import tsv_query_command


@pytest.mark.unit
def test_tsv_query_empty_filter(sample_tsv_file, caplog):
    """Test that empty filter returns all records with a warning."""
    with caplog.at_level(logging.WARNING):
        query_result, query_message = tsv_query_command(sample_tsv_file, filter="")

    assert len(query_result) == 5, "Empty filter should return all records"
    assert "Empty filter provided - returning ALL records" in caplog.text
    assert "ALL RECORDS" in query_message


@pytest.mark.unit
def test_tsv_query_none_filter(sample_tsv_file):
    """Test that None filter raises an appropriate error."""
    with pytest.raises(ValueError) as exc_info:
        tsv_query_command(sample_tsv_file, filter=None)

    assert "Filter cannot be None" in str(exc_info.value)


@pytest.mark.unit
def test_tsv_query_column_not_found(sample_tsv_file, caplog):
    """Test handling of filter with non-existent column."""
    with caplog.at_level(logging.WARNING):
        query_result, query_message = tsv_query_command(sample_tsv_file, filter="NonExistentColumn:Value")

    assert len(query_result) == 5, "Should return all records when column not found"
    assert "Column 'NonExistentColumn' not found in the TSV file" in caplog.text
    assert query_message == "Invalid filter conditions - returning ALL records"


@pytest.mark.unit
def test_tsv_query_value_not_found(sample_tsv_file, caplog):
    """Test handling of filter with non-existent values."""
    with caplog.at_level(logging.WARNING):
        query_result, query_message = tsv_query_command(sample_tsv_file, filter="Population:NonExistentPop")

    assert len(query_result) == 0, "Should return empty DataFrame when no values match"
    assert "None of the values ['NonExistentPop'] were found in column 'Population'" in caplog.text


@pytest.mark.unit
def test_tsv_query_invalid_filter_format(sample_tsv_file, caplog):
    """Test handling of incorrectly formatted filter."""
    with caplog.at_level(logging.WARNING):
        query_result, query_message = tsv_query_command(sample_tsv_file, filter="InvalidFormatNoColon")

    assert len(query_result) == 5, "Should return all records for invalid format"
    assert "Invalid filter format: 'InvalidFormatNoColon'" in caplog.text
    assert query_message == "Invalid filter conditions - returning ALL records"


@pytest.mark.unit
def test_tsv_query_filename_filter(sample_tsv_file):
    """Test filtering by Filename column."""
    query_result, _ = tsv_query_command(sample_tsv_file, filter="Filename:file1.vcf.gz")

    assert len(query_result) == 2, "Should return records with matching filename"
    assert list(query_result["Sample_ID"]) == ["S1", "S2"], "Should return correct samples"


@pytest.mark.unit
def test_tsv_query_multiple_conditions(sample_tsv_file):
    """Test filtering with multiple conditions (intersect query)."""
    query_result, _ = tsv_query_command(sample_tsv_file, filter="Population:Pop1;Sex:F")

    assert len(query_result) == 1, "Should return records matching both conditions"
    assert list(query_result["Sample_ID"]) == ["S2"], "Should return correct sample"


@pytest.mark.unit
def test_tsv_query_multiple_values(sample_tsv_file):
    """Test filtering with multiple values for a single column."""
    query_result, _ = tsv_query_command(sample_tsv_file, filter="Population:Pop1,Pop3")

    assert len(query_result) == 3, "Should return records matching any value"
    assert sorted(list(query_result["Sample_ID"])) == ["S1", "S2", "S5"], "Should return correct samples"


@pytest.mark.unit
def test_tsv_query_strip_hash_headers(tsv_with_hash_headers):
    """Test that column headers with # prefix are handled correctly."""
    query_result, _ = tsv_query_command(tsv_with_hash_headers, filter="Sample_ID:S1")

    assert len(query_result) == 1, "Should strip # from headers and match correctly"
    assert list(query_result.columns) == ["Sample_ID", "Population", "Filename"], "Should strip # from headers"
