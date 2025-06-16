"""
Tests for the "divbase-cli files" commands
"""

import shlex
from unittest.mock import patch

from typer.testing import CliRunner

from divbase_tools.divbase_cli import app

runner = CliRunner()


def test_list_files_command(mock_s3_manager, user_config_path, CONSTANTS):
    """Test the list_files command with a mocked S3FileManager."""
    command = f"files list --config {user_config_path}"

    with patch("divbase_tools.services.create_s3_file_manager", return_value=mock_s3_manager):
        result = runner.invoke(app, shlex.split(command))

    assert result.exit_code == 0

    for file in CONSTANTS["FILES_IN_BUCKET"]:
        assert file in result.stdout, f"File {file} not found in the output of the list_files command"
