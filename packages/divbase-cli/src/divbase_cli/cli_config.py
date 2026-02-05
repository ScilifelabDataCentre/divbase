"""
Settings for DivBase CLI.

This class creates a single 'settings' object at module load time that can be imported and used throughout the entire package.

The user config and tokens are stored in the user local app dir:
https://typer.tiangolo.com/tutorial/app-dir/
"""

import os
from dataclasses import dataclass
from pathlib import Path

import typer

APP_NAME = "divbase_cli"
APP_DIR = Path(typer.get_app_dir(APP_NAME))
CONFIG_PATH = APP_DIR / "config.yaml"
TOKENS_PATH = APP_DIR / ".secrets"
DEFAULT_METADATA_TSV_NAME = "sample_metadata.tsv"
DEFAULT_DIVBASE_API_URL = "http://localhost:8000/api"  # TODO - change to production URL when time comes


@dataclass
class DivBaseCLISettings:
    """
    Settings for DivBase CLI.

    NOTE: Do not create an instance of this class yourself,
    import the 'settings' instance created at this module's load time.
    TODO - TESTS could override the config and tokens path
    """

    CONFIG_PATH: Path = Path(os.getenv("DIVBASE_CONFIG_PATH", CONFIG_PATH))
    TOKENS_PATH: Path = Path(os.getenv("DIVBASE_TOKENS_PATH", TOKENS_PATH))
    DIVBASE_API_URL: str = os.getenv("DIVBASE_API_URL", DEFAULT_DIVBASE_API_URL)
    METADATA_TSV_NAME: str = os.getenv("DIVBASE_METADATA_TSV_NAME", DEFAULT_METADATA_TSV_NAME)
    LOGGING_ON: bool = bool(os.getenv("DIVBASE_LOGGING_ON", "True") == "True")
    LOG_LEVEL: str = os.getenv("DIVBASE_LOG_LEVEL", "INFO").upper()

    def __post_init__(self):
        valid_levels = ["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]
        if self.LOG_LEVEL not in valid_levels:
            raise ValueError(f"Invalid LOG_LEVEL: {self.LOG_LEVEL}. Must be one of {valid_levels}.")


cli_settings = DivBaseCLISettings()
