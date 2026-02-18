"""
Shared pytest fixtures for all unit tests.

Fixtures are defined here so they can be used across divbase_lib, divbase_cli,
and divbase_api unit tests without cross-file imports (which are unreliable under
pytest's --import-mode=importlib).
"""

import pytest


@pytest.fixture
def valid_tsv(tmp_path):
    """Simple valid TSV that passes all validation checks."""
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
def no_multi_values_tsv(tmp_path):
    """TSV with no semicolon-separated values in any cell."""
    tsv_content = """#Sample_ID\tPopulation\nS1\t1\nS2\t2\n"""
    tsv_file = tmp_path / "no_multi_values.tsv"
    tsv_file.write_text(tsv_content)
    return tsv_file


@pytest.fixture
def numeric_multi_values_tsv(tmp_path):
    """TSV with multi-value numeric cells and negative numbers."""
    tsv_content = """#Sample_ID\tScores\tValues\tTemperature\tLongitude\tLatitude\tElevation
S1\t1;2;3\t10;20\t-5.5\t-2.78305556\t51.5\t100
S2\t4;5\t30;40;50\t-10.2\t-0.12765\t52.2\t-50
S3\t6\t60\t0\t1.25\t50.8\t-100.5
S4\t7;8;9;10\t70\t15.5\t-3.5;-2.1\t49.5\t200
S5\t11\t80;90\t-20\t0\t48.2\t-25
"""
    tsv_file = tmp_path / "numeric_multi_values.tsv"
    tsv_file.write_text(tsv_content)
    return tsv_file


@pytest.fixture
def sample_tsv_with_numeric_data(tmp_path):
    """Comprehensive TSV with numeric, string, semicolon, float, and negative columns."""
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
def sample_tsv_with_mixed_type_column(tmp_path):
    """TSV with a mixed-type column (numeric-looking + non-numeric values)."""
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
    """TSV with a semicolon-separated column where one part is non-numeric ('1;1-2')."""
    tsv_content = """#Sample_ID\tCode\tPureNumericSemicolon\tWeight
S1\t1;1-2\t10;20;30\t12.5
S2\t3\t40\t18.8
S3\t5\t50;60\t15.0
S4\t7\t70;80;90\t20.0
"""
    tsv_file = tmp_path / "test_semicolon_mixed.tsv"
    tsv_file.write_text(tsv_content)
    return tsv_file


@pytest.fixture
def type_errors_tsv(tmp_path):
    """TSV with mixed-type columns and range/hyphen notation.

    Population: mixed int + string + range
    Test: mixed int + 'all'
    Code: pure string with number prefix
    Range: mixed int + hyphen-range notation
    """
    tsv_content = """#Sample_ID\tPopulation\tTest\tCode\tRange
S1\t1\t2\tA100\t1-2
S2\tabc\t3\tB200\t3
S3\t1;three;5\tall\tC300\t4
S4\t3-5\t4\tD400\t5
"""
    tsv_file = tmp_path / "type_errors.tsv"
    tsv_file.write_text(tsv_content)
    return tsv_file


@pytest.fixture
def array_notation_tsv(tmp_path):
    """TSV where one column uses '[1, 2, 3]' array notation instead of semicolons."""
    content = "#Sample_ID\tPopulation\tArea\n"
    content += "S1\t[1, 2, 3]\tNorth\n"
    content += "S2\t4\tEast\n"
    content += "S3\t5\tSouth\n"
    tsv_file = tmp_path / "array_notation.tsv"
    tsv_file.write_text(content)
    return tsv_file


@pytest.fixture
def array_notation_multiple_cols_tsv(tmp_path):
    """TSV where multiple columns use '[...]' array notation."""
    content = "#Sample_ID\tPopulation\tScores\n"
    content += "S1\t[1, 2]\t[10, 20, 30]\n"
    content += "S2\t3\t[40]\n"
    content += "S3\t5\t60\n"
    tsv_file = tmp_path / "array_notation_multi.tsv"
    tsv_file.write_text(content)
    return tsv_file


@pytest.fixture
def header_errors_tsv(tmp_path):
    """TSV with header errors: wrong first column, duplicate columns, empty column."""
    tsv_content = """SampleID\tPopulation\tArea\tArea\t
S1\t1\tNorth\tEast\tValue
"""
    tsv_file = tmp_path / "header_errors.tsv"
    tsv_file.write_text(tsv_content)
    return tsv_file


@pytest.fixture
def sample_errors_tsv(tmp_path):
    """TSV with Sample_ID errors: empty, semicolons, duplicates."""
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
    """TSV with formatting errors: wrong column count, commas, whitespace."""
    tsv_content = """#Sample_ID\tPopulation\tArea
S1\t1\tNorth
S2\t2,3\tEast
S3 \t 4 \t West 
S4\t5
"""
    tsv_file = tmp_path / "format_errors.tsv"
    tsv_file.write_text(tsv_content)
    return tsv_file


@pytest.fixture
def sample_tsv_with_invalid_sample_ids(tmp_path):
    """TSV with empty and duplicate Sample_IDs (both raise errors during load)."""
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
    """TSV that omits the Sample_ID column entirely."""
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
    """TSV with duplicate Sample_IDs (raises error during load)."""
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
def sample_tsv_with_edge_cases(tmp_path):
    """TSV with edge cases: unicode, hyphens in strings, and whitespace in Sample_IDs.

    NOTE: S2 and S3 have leading/trailing whitespace in Sample_ID — the validator
    strips these, so the exported file will DIFFER from the original. Do NOT use
    this fixture for identity roundtrip tests.
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
def negative_numeric_columns():
    """Column names in numeric_multi_values_tsv that contain negative numbers."""
    return ["Temperature", "Longitude", "Latitude", "Elevation"]
