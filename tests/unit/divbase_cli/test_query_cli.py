"""Unit tests for query CLI sample input normalization."""

import pytest
import typer

from divbase_cli.cli_commands.query_cli import _normalize_samples_input


def test_normalize_samples_option_splits_comma_separated_values():
    "Test that the --samples option correctly splits comma-separated sample IDs and strips whitespace."
    normalized, warnings = _normalize_samples_input(samples="S1, S2 ,S3", samples_file=None)
    assert normalized == ["S1", "S2", "S3"]
    assert warnings == []


def test_normalize_samples_file_ignores_comments_and_blank_lines(tmp_path):
    """Test that the --samples-file option correctly ignores blank lines and lines starting with #, and normalizes the remaining sample IDs."""
    samples_file = tmp_path / "samples.txt"
    samples_file.write_text(
        "\n".join(
            [
                "# Sample IDs for run",
                "",
                "S1",
                "  # temporary note",
                "  S2  ",
            ]
        ),
        encoding="utf-8",
    )

    normalized, warnings = _normalize_samples_input(samples=None, samples_file=samples_file)
    assert normalized == ["S1", "S2"]
    assert warnings == []


def test_normalize_samples_file_warns_when_only_one_sample_is_present(tmp_path):
    """Test that the --samples-file option issues a warning when only one sample ID is found after ignoring blank/comment lines, as this may indicate a formatting issue with the file."""
    samples_file = tmp_path / "single_sample.txt"
    samples_file.write_text("# single sample run\nS1\n", encoding="utf-8")

    normalized, warnings = _normalize_samples_input(samples=None, samples_file=samples_file)
    assert normalized == ["S1"]
    assert any("Only one sample ID was found" in warning for warning in warnings)


def test_normalize_samples_file_errors_for_delimited_lines(tmp_path):
    """Test that the --samples-file option raises an error if any line contains delimiters, as this indicates an invalid format (multiple sample IDs on one line instead of one per line)."""
    samples_file = tmp_path / "suspicious_samples.txt"
    samples_file.write_text(
        "\n".join(
            [
                "S1,S2,S3",
                "S4;S5",
                "S6|S7",
                "S8\tS9",
            ]
        ),
        encoding="utf-8",
    )

    with pytest.raises(typer.BadParameter, match="expected one sample ID per line with no delimiters"):
        _normalize_samples_input(samples=None, samples_file=samples_file)


def test_normalize_samples_file_error_lists_found_delimiters(tmp_path):
    """Test that the error message from --samples-file option includes which delimiters were found and on which lines, to help users identify formatting issues with their samples file."""
    samples_file = tmp_path / "suspicious_delimiters.txt"
    samples_file.write_text("S1;S2\nS3|S4\n", encoding="utf-8")

    with pytest.raises(typer.BadParameter) as excinfo:
        _normalize_samples_input(samples=None, samples_file=samples_file)

    error_message = str(excinfo.value)
    assert "';'" in error_message
    assert "'|'" in error_message


def test_normalize_samples_file_errors_when_only_comments_and_blank_lines(tmp_path):
    """Test that the --samples-file option raises an error when the file contains only comments and blank lines, as this indicates an invalid format."""
    samples_file = tmp_path / "empty_samples.txt"
    samples_file.write_text("# no data yet\n\n  # still no data\n", encoding="utf-8")

    with pytest.raises(typer.BadParameter, match="has no sample IDs"):
        _normalize_samples_input(samples=None, samples_file=samples_file)
