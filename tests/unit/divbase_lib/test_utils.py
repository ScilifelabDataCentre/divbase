"""Unit test for the utils module in the divbase_cli package."""

from datetime import datetime, timezone
from zoneinfo import ZoneInfo

import pytest

from divbase_lib.utils import (
    format_datetime,
    format_file_size,
    split_semicolon_bcftools_command_segments,
    to_unix_timestamp,
)


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


@pytest.mark.parametrize(
    "dt, expected_unix_timestamp",
    [
        (datetime(2026, 1, 1, 0, 0, 0, tzinfo=timezone.utc), 1767225600),
        (datetime(2026, 1, 1, 1, 0, 0, tzinfo=ZoneInfo("Europe/Stockholm")), 1767225600),
        (datetime(2026, 1, 1, 0, 0, 0), 1767225600),  # naive datetimes should be assumed to be UTC
    ],
)
def test_to_unix_timestamp(dt, expected_unix_timestamp):
    """Test that to_unix_timestamp converts a datetime to the correct unix timestamp."""
    assert to_unix_timestamp(dt) == expected_unix_timestamp


def test_format_datetime():
    """Test that format_datetime displays Europe/Stockholm time by default."""
    dt = datetime(2026, 1, 1, 0, 0, 0, tzinfo=timezone.utc)
    naive_dt = datetime(2026, 1, 1, 0, 0, 0)  # should be treated as UTC if not specified.
    expected_dt = "2026-01-01 01:00:00 CET"

    assert format_datetime(dt) == expected_dt
    assert format_datetime(naive_dt) == expected_dt


# NOTE, we don't test format_datetime_for_cli or format_unix_timestamp_for_cli as they are
# simple wrappers and we would need to handle local timezones diffs in the test, so not worth the hassle.


@pytest.mark.parametrize(
    "command,expected_segments",
    [
        ("view -s; view -r 1", ["view -s", " view -r 1"]),
        ("view -i 'FILTER=\"A;B\"'; view -r 1", ["view -i 'FILTER=\"A;B\"'", " view -r 1"]),
        ('view -i "FILTER=\\"A;B\\""; view -r 1', ['view -i "FILTER=\\"A;B\\""', " view -r 1"]),
        ("view -r 1;", ["view -r 1", ""]),
    ],
)
def test_split_semicolon_command_segments(command, expected_segments):
    """Test that command splitting keeps semicolons inside quoted substrings."""
    assert split_semicolon_bcftools_command_segments(command) == expected_segments
