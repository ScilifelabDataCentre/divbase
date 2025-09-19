"""
Config for DivBase API.

This module creates a single instance of the Settings class,
which other modules can import directly to access configuration settings.

For local development these settings are defined by environment variables,
which are set in the docker compose file.
"""

import os
from dataclasses import dataclass, field

from pydantic import SecretStr


@dataclass
class DBSettings:
    """PostgreSQL database configuration settings."""

    url: SecretStr = SecretStr(os.getenv("POSTGRES_URL", "NOT_SET"))
    echo_db_output: bool = bool(os.getenv("DEBUG", "False") == "True")  # anything but "True" is considered False


@dataclass
class FlowerSettings:
    """Flower configuration settings."""

    user: str = os.getenv("FLOWER_USER", "NOT_SET")
    password: SecretStr = SecretStr(os.getenv("FLOWER_PASSWORD", "NOT_SET"))
    url: str = os.getenv("FLOWER_URL", "http://flower:5555")


@dataclass
class Settings:
    """Configuration settings for DivBase API."""

    database: DBSettings = field(default_factory=DBSettings)
    flower: FlowerSettings = field(default_factory=FlowerSettings)

    def __post_init__(self):
        """
        Validate all required settings are actually set.
        This means that later on in the codebase we don't have to check for any non set values, we can just assume they are set.
        """
        required_fields = [self.flower.url, self.flower.user]
        for setting in required_fields:
            if setting == "NOT_SET":
                raise ValueError(f"A required environment variable was not set: {setting}")

        required_secret_fields = [self.database.url, self.flower.password]
        for secret in required_secret_fields:
            if secret.get_secret_value() == "NOT_SET":
                # TODO - check if outputted error show the name of the missing variable or just *****
                raise ValueError(f"A required secret environment variable was not set: {secret}")


# This instance can be imported and used throughout the application.
settings = Settings()
