"""
Unit tests for validating CLI versions when sending requests to the API.
"""

import pytest

from divbase_api.services.validate_cli_versions import cli_update_available, cli_version_outdated

outdated_cases = [
    # (cli_version, minimum_cli_version, expected_result)
    ("3.1.0", "3.1.0", False),
    ("3.1.0", "3.2.0", True),
    ("3.1.0", "4.0.0", True),
    ("3.1.0", "2.4.9", False),
    ("3.1.0", "2.5.0", False),
    ("3.1.0", "2.4.0", False),
    ("3.1.0", "2.9.9", False),
    ("3.1.0", "3.0.0", False),
    ("3.1.0", "3.0.9", False),
    ("3.1.0", "3.1.1", True),
    ("3.1.0", "3.1.0rc1", False),
    ("2.1.0rc123", "2.1.0", False),
    ("2.1.0", "2.1.0rc123", False),
    ("2.1.0", "2.1.1", True),
    ("2.1.1", "2.1.0", False),
    ("1.0.0", "1.0.0", False),
    ("1.0.0", "0.9.9", False),
    ("1.0.0", "1.0.1", True),
    ("1.0.0", "1.1.0", True),
    ("1.0.0", "2.0.0", True),
]


@pytest.mark.parametrize("cli_version, minimum_cli_version, expected_result", outdated_cases)
def test_cli_version_outdated(cli_version: str, minimum_cli_version: str, expected_result: bool):
    result = cli_version_outdated(cli_version=cli_version, minimum_cli_version=minimum_cli_version)
    assert result == expected_result


update_cases = [
    # (cli_version, latest_cli_version, expected_result)
    ("3.1.0", "3.1.0", False),
    ("3.1.0", "3.2.0", True),
    ("3.1.0", "4.0.0", True),
    ("3.1.0", "2.4.9", False),
    ("3.1.0", "2.5.0", False),
    ("3.1.0", "2.4.0", False),
    ("3.1.0", "2.9.9", False),
    ("3.1.0", "3.0.0", False),
    ("3.1.0", "3.0.9", False),
    ("3.1.0", "3.1.1", True),
    ("3.1.0", "3.1.0rc1", False),
    ("2.1.0rc123", "2.1.0", False),
    ("2.1.0", "2.1.0rc123", False),
    ("2.1.0", "2.1.1", True),
    ("2.1.1", "2.1.0", False),
    ("1.0.0", "1.0.0", False),
    ("1.0.0", "0.9.9", False),
    ("1.0.0", "1.0.1", True),
    ("1.0.0", "1.1.0", True),
    ("1.0.0", "2.0.0", True),
]


@pytest.mark.parametrize("cli_version, latest_cli_version, expected_result", update_cases)
def test_cli_update_available(cli_version: str, latest_cli_version: str, expected_result: bool):
    result = cli_update_available(cli_version=cli_version, latest_cli_version=latest_cli_version)
    assert result == expected_result
