import logging
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

from divbase_tools.exceptions import SidecarColumnNotFoundError, SidecarInvalidFilterError, SidecarNoDataLoadedError
from divbase_tools.queries import SidecarQueryManager


@pytest.mark.unit
def test_sidecar_manager_SidecarNoDataLoadedError_for_non_existent_file(create_sidecar_manager):
    """
    Test that SidecarQueryManager raises SidecarNoDataLoadedError when the file path does not exist.
    """
    with pytest.raises(SidecarNoDataLoadedError):
        create_sidecar_manager("non_existent_file.tsv")


@pytest.mark.unit
@patch("divbase_tools.queries.SidecarQueryManager.load_file")
def test_sidecar_manager_no_data_loaded(mock_load_file, create_sidecar_manager):
    """
    Test that SidecarNoDataLoadedError is raised when no data (self.df=None) is loaded in the SidecarQueryManager.
    """
    mock_load_file.return_value = None
    sidecar_manager = create_sidecar_manager("sample_file.tsv")

    with pytest.raises(SidecarNoDataLoadedError):
        sidecar_manager.run_query()


@pytest.mark.unit
def test_run_query_invalid_filter_raises_custom_exception(create_sidecar_manager):
    """Test that SidecarQueryManager raises SidecarInvalidFilterError
    when an invalid filter format is provided."""

    with tempfile.NamedTemporaryFile("w", delete=False, suffix=".tsv") as tmp_valid_tsv:
        tmp_valid_tsv.write("col1\tcol2\nA\t1\nB\t2\n")
        tmp_path = Path(tmp_valid_tsv.name)

    manager = create_sidecar_manager(tmp_path)
    invalid_filter = "col1A,B"

    with pytest.raises(SidecarInvalidFilterError) as excinfo:
        manager.run_query(invalid_filter)

    assert "Invalid filter format" in str(excinfo.value)
    assert invalid_filter in str(excinfo.value)

    tmp_path.unlink()


@pytest.mark.unit
def test_empty_string(valid_tsv_path, caplog):
    """
    Test that an empty input filter returns all records and that unique values are extracted correctly.
    """
    manager = SidecarQueryManager(valid_tsv_path)
    with caplog.at_level("WARNING"):
        manager.run_query("")

    assert "Empty filter provided - returning ALL records. This may be a large result set." in caplog.text
    assert sorted(manager.get_unique_values("col1")) == ["A", "B"]
    assert sorted(manager.get_unique_values("col2")) == [1, 2, 3]


@pytest.mark.unit
def test_get_unique_values_no_query_result(valid_tsv_path):
    """
    Test that get_unique_values raises SidecarColumnNotFoundError when no query result is available.
    To achieve this, the test does not call run_query before calling get_unique_values (i.e. query_result=None).
    """
    manager = SidecarQueryManager(valid_tsv_path)

    with pytest.raises(SidecarColumnNotFoundError) as excinfo:
        manager.get_unique_values("col1")
    assert "No query result available" in str(excinfo.value)


@pytest.mark.unit
def test_get_unique_values_column_not_found(valid_tsv_path):
    """
    Test that get_unique_values raises SidecarColumnNotFoundError if a non-existent column is requested.
    Do this by running a valid query first, then trying to get unique values for a column that does not exist.
    """
    manager = SidecarQueryManager(valid_tsv_path)
    manager.run_query("col1:A,C")

    with pytest.raises(SidecarColumnNotFoundError) as excinfo:
        manager.get_unique_values("nonexistent_column")
    assert "not found in query result" in str(excinfo.value)


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
def test_tsv_query_empty_filter(sample_tsv_file, caplog, create_sidecar_manager):
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
def test_tsv_query_none_filter(sample_tsv_file, create_sidecar_manager):
    """Test that None filter raises an appropriate error."""
    manager = create_sidecar_manager(sample_tsv_file)
    with pytest.raises(SidecarInvalidFilterError) as exc_info:
        manager.run_query(filter_string=None)

    assert "Filter cannot be None" in str(exc_info.value)


@pytest.mark.unit
def test_tsv_query_column_not_found(sample_tsv_file, caplog, create_sidecar_manager):
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
def test_tsv_query_value_not_found(sample_tsv_file, caplog, create_sidecar_manager):
    """Test handling of filter with non-existent values."""
    with caplog.at_level(logging.WARNING):
        manager = create_sidecar_manager(sample_tsv_file)
        manager.run_query(filter_string="Population:NonExistentPop")
        query_result = manager.query_result

    assert len(query_result) == 0, "Should return empty DataFrame when no values match"
    assert "None of the values ['NonExistentPop'] were found in column 'Population'" in caplog.text


@pytest.mark.unit
def test_tsv_query_invalid_filter_format(sample_tsv_file, caplog, create_sidecar_manager):
    """Test handling of incorrectly formatted filter."""
    with caplog.at_level(logging.WARNING):
        manager = create_sidecar_manager(sample_tsv_file)
        with pytest.raises(SidecarInvalidFilterError) as exc_info:
            manager.run_query(filter_string="InvalidFormatNoColon")
        assert "Invalid filter format" in str(exc_info.value)


@pytest.mark.unit
def test_tsv_query_filename_filter(sample_tsv_file, create_sidecar_manager):
    """Test filtering by Filename column."""
    manager = create_sidecar_manager(sample_tsv_file)
    manager.run_query(filter_string="Filename:file1.vcf.gz")
    query_result = manager.query_result

    assert len(query_result) == 2, "Should return records with matching filename"
    assert list(query_result["Sample_ID"]) == ["S1", "S2"], "Should return correct samples"


@pytest.mark.unit
def test_tsv_query_multiple_conditions(sample_tsv_file, create_sidecar_manager):
    """Test filtering with multiple conditions (intersect query)."""
    manager = create_sidecar_manager(sample_tsv_file)
    manager.run_query(filter_string="Population:Pop1;Sex:F")
    query_result = manager.query_result

    assert len(query_result) == 1, "Should return records matching both conditions"
    assert list(query_result["Sample_ID"]) == ["S2"], "Should return correct sample"


@pytest.mark.unit
def test_tsv_query_multiple_values(sample_tsv_file, create_sidecar_manager):
    """Test filtering with multiple values for a single column."""
    manager = create_sidecar_manager(sample_tsv_file)
    manager.run_query(filter_string="Population:Pop1,Pop3")
    query_result = manager.query_result

    assert len(query_result) == 3, "Should return records matching any value"
    assert sorted(list(query_result["Sample_ID"])) == ["S1", "S2", "S5"], "Should return correct samples"


@pytest.mark.unit
def test_tsv_query_strip_hash_headers(tsv_with_hash_headers, create_sidecar_manager):
    """Test that column headers with # prefix are handled correctly."""
    manager = create_sidecar_manager(tsv_with_hash_headers)
    manager.run_query(filter_string="Sample_ID:S1")
    query_result = manager.query_result

    assert len(query_result) == 1, "Should strip # from headers and match correctly"
    assert list(query_result.columns) == ["Sample_ID", "Population", "Filename"], "Should strip # from headers"
