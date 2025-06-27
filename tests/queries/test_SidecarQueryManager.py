import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest
from conftest import create_sidecar_manager

from divbase_tools.exceptions import SidecarColumnNotFoundError, SidecarInvalidFilterError, SidecarNoDataLoadedError
from divbase_tools.queries import SidecarQueryManager


@pytest.mark.unit
def test_sidecar_manager_SidecarNoDataLoadedError_for_non_existent_file():
    """
    Test that SidecarQueryManager raises SidecarNoDataLoadedError when the file path does not exist.
    """
    with pytest.raises(SidecarNoDataLoadedError):
        create_sidecar_manager("non_existent_file.tsv")


@pytest.mark.unit
@patch("divbase_tools.queries.SidecarQueryManager.load_file")
def test_sidecar_manager_no_data_loaded(mock_load_file):
    """
    Test that SidecarNoDataLoadedError is raised when no data (self.df=None) is loaded in the SidecarQueryManager.
    """
    mock_load_file.return_value = None
    sidecar_manager = create_sidecar_manager("sample_file.tsv")

    with pytest.raises(SidecarNoDataLoadedError):
        sidecar_manager.run_query()


@pytest.mark.unit
def test_run_query_invalid_filter_raises_custom_exception():
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
