"""
Settings for DivBase CLI.

This class creates a single 'settings' object at module load time that can be imported and used throughout the entire package.
"""

import os
from dataclasses import dataclass
from pathlib import Path

from dotenv import load_dotenv

DEFAULT_CONFIG_PATH = Path.home() / ".config" / "divbase" / "config.yaml"
DEFAULT_TOKENS_PATH = Path.home() / ".config" / "divbase" / ".secrets"
DEFAULT_METADATA_TSV_NAME = "sample_metadata.tsv"
DEFAULT_DIVBASE_API_URL = "http://localhost:8000/api"  # TODO - change to production URL when time comes
DEFAULT_S3_URL = "http://localhost:9000"  # TODO - change to production URL when time comes


# This is used to get pytest to read the enviroment variables before the cli_settings instance is created.
# Pytest needs to overwrite the default paths and urls to prevent test pollution.
# Trying to set these env variables in conftest.py does not work
# as this module is imported before conftest.py is processed.
if os.getenv("PYTEST", "").lower() in ("1", "true"):
    load_dotenv("tests/cli_settings.txt")


@dataclass
class DivBaseCLISettings:
    """
    Settings for DivBase CLI.
    NOTE: Do not create an instance of this class yourself,
    import the 'settings' instance created at this module's load time.

    For end user and local development, the defaults are likely fine, so no need to set these environment variables.
    The use case for changing the env variables is for running with pytest.
    """

    CONFIG_PATH: Path = Path(os.getenv("DIVBASE_CONFIG_PATH", DEFAULT_CONFIG_PATH))
    TOKENS_PATH: Path = Path(os.getenv("DIVBASE_TOKENS_PATH", DEFAULT_TOKENS_PATH))
    DIVBASE_API_URL: str = os.getenv("DIVBASE_API_URL", DEFAULT_DIVBASE_API_URL)
    S3_URL: str = os.getenv("DIVBASE_S3_URL", DEFAULT_S3_URL)
    METADATA_TSV_NAME: str = os.getenv("DIVBASE_METADATA_TSV_NAME", DEFAULT_METADATA_TSV_NAME)


cli_settings = DivBaseCLISettings()
