"""
Tests for the "divbase-cli files" commands

All tests are run against a MinIO server on localhost from docker-compose.
"""

import shlex
from pathlib import Path

from typer.testing import CliRunner

from divbase_tools.divbase_cli import app

runner = CliRunner()

FIXTURES_DIR = Path(__file__).parent.parent / "fixtures"


def test_list_files_command(user_config_path, CONSTANTS):
    """Test the list_files command"""
    command = f"files list --config {user_config_path}"

    result = runner.invoke(app, shlex.split(command))
    assert result.exit_code == 0

    for file in CONSTANTS["FILES_IN_BUCKET"]:
        assert file in result.stdout, f"File {file} not found in the output of the list_files command"
