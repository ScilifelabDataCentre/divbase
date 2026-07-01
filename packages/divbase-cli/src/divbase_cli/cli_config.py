"""
Settings for DivBase CLI.

This class creates a single 'settings' object at module load time that can be imported and used throughout the entire package.

The user config and tokens are stored in the users local app dir:
https://typer.tiangolo.com/tutorial/app-dir/
"""

import os
from dataclasses import dataclass, field
from pathlib import Path

import typer
from pydantic import SecretStr

APP_NAME = "divbase-cli"
APP_DIR = Path(typer.get_app_dir(APP_NAME))
DEFAULT_METADATA_TSV_NAME = "sample_metadata.tsv"
DEFAULT_LOG_LEVEL = "INFO"
DEV_MODE = os.getenv("DIVBASE_DEV", "0") == "1"
DEFAULT_CONFIG_PATH = APP_DIR / "config.yaml"
DEFAULT_TOKENS_PATH = APP_DIR / ".secrets"
DEFAULT_PATS_PATH = APP_DIR / ".pat"  # no plural because only 1 stored.

if DEV_MODE:
    DEFAULT_DIVBASE_API_URL = "http://localhost:8000/api"
    DEFAULT_LOGGING_ON = "1"
    DEFAULT_TRACEBACKS_ON = "1"
else:
    DEFAULT_DIVBASE_API_URL = "https://divbase.scilifelab-2-prod.sys.kth.se/api"
    DEFAULT_LOGGING_ON = "0"
    DEFAULT_TRACEBACKS_ON = "0"


@dataclass
class DivBaseCLISettings:
    """
    Settings for DivBase CLI.

    You do not need to create an instance of this class yourself,
    instead, import the 'cli_settings' instance created at this module's load time.
    """

    CONFIG_PATH: Path = Path(os.getenv("DIVBASE_CLI_CONFIG_PATH", DEFAULT_CONFIG_PATH))

    # for tokens stored via OS keyring, we use these to define a unique lookup key.
    KEYRING_SERVICE: str = os.getenv("DIVBASE_KEYRING_SERVICE", "divbase-cli")
    KEYRING_TOKENS_USERNAME: str = "tokens"
    KEYRING_PATS_USERNAME: str = "pat"
    # Fallback paths for tokens (JWTs) and PATs storage if keyring cannot be used on the device.
    TOKENS_FALLBACK_PATH: Path = Path(os.getenv("DIVBASE_CLI_TOKENS_PATH", DEFAULT_TOKENS_PATH))
    PATS_FALLBACK_PATH: Path = Path(os.getenv("DIVBASE_CLI_PATS_PATH", DEFAULT_PATS_PATH))

    DIVBASE_API_URL: str = os.getenv("DIVBASE_API_URL", DEFAULT_DIVBASE_API_URL)
    METADATA_TSV_NAME: str = os.getenv("DIVBASE_METADATA_TSV_NAME", DEFAULT_METADATA_TSV_NAME)
    TRACEBACKS_ON: bool = os.getenv("DIVBASE_TRACEBACKS_ON", DEFAULT_TRACEBACKS_ON) == "1"
    LOGGING_ON: bool = os.getenv("DIVBASE_LOGGING_ON", DEFAULT_LOGGING_ON) == "1"
    LOG_LEVEL: str = os.getenv("DIVBASE_LOG_LEVEL", DEFAULT_LOG_LEVEL).upper()

    DIVBASE_API_PAT: SecretStr | None = field(init=False)

    def __post_init__(self):
        if os.getenv("DIVBASE_API_PAT"):
            self.DIVBASE_API_PAT = SecretStr(os.environ["DIVBASE_API_PAT"])
        else:
            self.DIVBASE_API_PAT = None

        valid_levels = ["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]
        if self.LOG_LEVEL not in valid_levels:
            raise ValueError(f"Invalid LOG_LEVEL: {self.LOG_LEVEL}. Must be one of {valid_levels}.")

        if self.DIVBASE_API_URL.endswith("/"):
            self.DIVBASE_API_URL = self.DIVBASE_API_URL[:-1]


cli_settings = DivBaseCLISettings()
