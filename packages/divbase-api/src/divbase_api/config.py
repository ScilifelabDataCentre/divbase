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
class APISettings:
    """API configuration settings."""

    log_level: str = os.getenv("LOG_LEVEL", "INFO").upper()
    first_admin_email: str = os.getenv("FIRST_ADMIN_EMAIL", "NOT_SET")
    first_admin_password: SecretStr = SecretStr(os.getenv("FIRST_ADMIN_PASSWORD", "NOT_SET"))
    worker_service_email: str = os.getenv("WORKER_SERVICE_EMAIL", "NOT_SET")
    worker_service_password: SecretStr = SecretStr(os.getenv("WORKER_SERVICE_PASSWORD", "NOT_SET"))


@dataclass
class DBSettings:
    """PostgreSQL database configuration settings."""

    url: SecretStr = SecretStr(os.getenv("DATABASE_URL", "NOT_SET"))
    echo_db_output: bool = bool(os.getenv("DB_ECHO", "False") == "True")  # anything but "True" is considered False


@dataclass
class FlowerSettings:
    """Flower configuration settings."""

    user: str = os.getenv("FLOWER_USER", "NOT_SET")
    password: SecretStr = SecretStr(os.getenv("FLOWER_PASSWORD", "NOT_SET"))
    url: str = os.getenv("FLOWER_URL", "http://flower:5555")


@dataclass
class S3Settings:
    """
    S3 configuration settings.

    External is used by the API to generate pre-signed URLs for end users.
    Internal is used by the backend to communicate with S3 directly.
    (This wont be needed later when S3 is moved outside of local docker network).
    """

    s3_external_url: str = os.getenv("S3_EXTERNAL_URL", "http://localhost:9000")
    s3_internal_url: str = os.getenv("S3_INTERNAL_URL", "http://minio:9000")
    access_key: SecretStr = SecretStr(os.getenv("S3_ACCESS_KEY", "NOT_SET"))
    secret_key: SecretStr = SecretStr(os.getenv("S3_SECRET_KEY", "NOT_SET"))


@dataclass
class JWTSettings:
    """JSON Web Token (JWT) configuration settings."""

    secret_key: SecretStr = SecretStr(os.getenv("JWT_SECRET_KEY", "NOT_SET"))
    algorithm: str = os.getenv("JWT_ALGORITHM", "HS256")
    access_token_expires_seconds: int = int(os.getenv("JWT_ACCESS_EXPIRES_SECONDS", 15 * 60))  # 15 mins
    refresh_token_expires_seconds: int = int(os.getenv("JWT_REFRESH_EXPIRES_SECONDS", 60 * 60 * 24 * 7))  # 7 days


@dataclass
class Settings:
    """Configuration settings for DivBase API."""

    api: APISettings = field(default_factory=APISettings)
    database: DBSettings = field(default_factory=DBSettings)
    flower: FlowerSettings = field(default_factory=FlowerSettings)
    s3: S3Settings = field(default_factory=S3Settings)
    jwt: JWTSettings = field(default_factory=JWTSettings)

    def __post_init__(self):
        """
        Validate all required settings are actually set.
        This means that later on in the codebase we don't have to check for any non set values, we can just assume they are set.
        """
        required_fields = {
            "DATABASE_URL": self.database.url,
            "FLOWER_URL": self.flower.url,
            "FLOWER_USER": self.flower.user,
            "FLOWER_PASSWORD": self.flower.password,
            "JWT_SECRET_KEY": self.jwt.secret_key,
            "S3_ACCESS_KEY": self.s3.access_key,
            "S3_SECRET_KEY": self.s3.secret_key,
        }
        for setting_name, setting in required_fields.items():
            if isinstance(setting, SecretStr) and setting.get_secret_value() == "NOT_SET":
                raise ValueError(f"A required secret environment variable was not set: {setting_name}")
            if isinstance(setting, str) and setting == "NOT_SET":
                raise ValueError(f"A required environment variable was not set: {setting_name}")


# This instance can be imported and used throughout the application.
settings = Settings()
