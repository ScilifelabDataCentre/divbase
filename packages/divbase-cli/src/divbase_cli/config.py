"""
Settings for DivBase CLI.

This class creates a single 'settings' object at module load time that can be imported and used throughout the entire package.
"""

from dataclasses import dataclass
from pathlib import Path


@dataclass
class DivBaseCLISettings:
    """
    Settings for DivBase CLI.
    """

    DEFAULT_CONFIG_PATH: Path = Path.home() / ".config" / "divbase" / "config.yaml"
    DEFAULT_TOKEN_PATH: Path = Path.home() / ".config" / "divbase" / ".env"
    DEFAULT_DIVBASE_API_URL: str = "http://localhost:8000/api/v1"  # TODO - change to production URL when time comes
    DEFAULT_S3_URL: str = "http://localhost:9000"  # TODO - change to production URL when time comes
    DEFAULT_METADATA_TSV_NAME: str = "sample_metadata.tsv"


settings = DivBaseCLISettings()
