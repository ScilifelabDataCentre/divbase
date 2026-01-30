"""Unit test for the utils module in the divbase_cli package."""

import pytest

from divbase_cli.utils import format_file_size


@pytest.mark.parametrize(
    "size_bytes, expected_output",
    [
        (None, "N/A"),
        (0, "0 B"),
        (500, "500.00 B"),
        (999, "999.00 B"),
        (1000, "1.00 KB"),
        (1500, "1.50 KB"),
        (999_999, "1000.00 KB"),
        (1_000_000, "1.00 MB"),
        (2_500_000, "2.50 MB"),
        (1_000_000_000, "1.00 GB"),
        (3_800_000_000, "3.80 GB"),
        (1_000_000_000_000, "1.00 TB"),
        (4_200_000_000_000, "4.20 TB"),
        (1_234_567_890_123, "1.23 TB"),
        (1023, "1.02 KB"),
        (1000.5, "1.00 KB"),
    ],
)
def test_format_file_size(size_bytes, expected_output):
    """
    Test that format_file_size correctly converts byte sizes to human-readable strings.
    """
    assert format_file_size(size_bytes) == expected_output
