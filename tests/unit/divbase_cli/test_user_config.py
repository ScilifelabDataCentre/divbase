from pathlib import Path

import pytest

from divbase_cli.user_config import load_user_config


def test_show_user_config_fails_if_no_config_file(tmp_path: Path) -> None:
    """Test that the show command fails if no config file exists."""

    with pytest.raises(FileNotFoundError):
        load_user_config(config_path=tmp_path / "non_existent_config.yaml")
