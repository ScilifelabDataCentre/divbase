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

    environment: str = os.getenv("DIVBASE_ENV", "NOT_SET")
    frontend_base_url: str = os.getenv("FRONTEND_BASE_URL", "NOT_SET")
    log_level: str = os.getenv("LOG_LEVEL", "INFO").upper()
    first_admin_email: str = os.getenv("FIRST_ADMIN_EMAIL", "NOT_SET")
    first_admin_password: SecretStr = SecretStr(os.getenv("FIRST_ADMIN_PASSWORD", "NOT_SET"))


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

    - S3_endpoint_url is for DivBase Server to directly interact with S3.
    - S3_presigning_url is used to generate pre-signed URLs for end users.
    These may or may not be the same URL depending on the setup.
    """

    endpoint_url: str = os.getenv("S3_ENDPOINT_URL", "NOT_SET")
    presigning_url: str = os.getenv("S3_PRESIGNING_URL", "NOT_SET")
    access_key: SecretStr = SecretStr(os.getenv("S3_SERVICE_ACCOUNT_ACCESS_KEY", "NOT_SET"))
    secret_key: SecretStr = SecretStr(os.getenv("S3_SERVICE_ACCOUNT_SECRET_KEY", "NOT_SET"))


@dataclass
class JWTSettings:
    """JSON Web Token (JWT) configuration settings."""

    secret_key: SecretStr = SecretStr(os.getenv("JWT_SECRET_KEY", "NOT_SET"))
    algorithm: str = os.getenv("JWT_ALGORITHM", "HS256")
    access_token_expires_seconds: int = int(os.getenv("JWT_ACCESS_EXPIRES_SECONDS", 15 * 60))  # 15 mins
    refresh_token_expires_seconds: int = int(os.getenv("JWT_REFRESH_EXPIRES_SECONDS", 60 * 60 * 24 * 7))  # 7 days


@dataclass
class EmailSettings:
    """
    Email configuration settings.
    Currently only working for local development with Mailpit
    """

    smtp_server: str = field(init=False)
    smtp_port: int = field(init=False)
    smtp_tls: bool = field(init=False)
    smtp_ssl: bool = field(init=False)

    smtp_user: str | None = os.getenv("SMTP_USER", None)
    smtp_password: SecretStr | None = field(init=False)

    from_email: str = os.getenv("FROM_EMAIL", "noreply-divbase@scilifelab.se")

    # token expiration times included in emails
    email_verify_expires_seconds: int = int(os.getenv("EMAIL_VERIFY_EXPIRES_SECONDS", 60 * 60 * 24))  # 24 hours
    password_reset_expires_seconds: int = int(os.getenv("PASSWORD_RESET_EXPIRES_SECONDS", 60 * 60))  # 1 hour

    def __post_init__(self):
        """Handle enviroment specific email settings."""
        if os.getenv("DIVBASE_ENV") in ["local_dev", "test"]:
            # using mailpit in docker stack
            self.smtp_server = "mailpit"
            self.smtp_port = 1025
            self.smtp_password = None
            self.smtp_tls = False
            self.smtp_ssl = False

        else:
            self.smtp_server = os.getenv("SMTP_SERVER", "smtp-relay.gmail.com")
            self.smtp_port = int(os.getenv("SMTP_PORT", 587))

            self.smtp_tls = bool(os.getenv("SMTP_TLS", "True") == "True")
            self.smtp_ssl = bool(os.getenv("SMTP_SSL", "False") == "True")
            if self.smtp_tls and self.smtp_ssl:
                raise ValueError("SMTP_TLS and SMTP_SSL cannot both be True.")

            smtp_password = os.getenv("SMTP_PASSWORD", None)
            if smtp_password:
                self.smtp_password = SecretStr(smtp_password)
            else:
                self.smtp_password = None


@dataclass
class Settings:
    """Configuration settings for DivBase API."""

    api: APISettings = field(default_factory=APISettings)
    database: DBSettings = field(default_factory=DBSettings)
    flower: FlowerSettings = field(default_factory=FlowerSettings)
    s3: S3Settings = field(default_factory=S3Settings)
    jwt: JWTSettings = field(default_factory=JWTSettings)
    email: EmailSettings = field(default_factory=EmailSettings)

    def validate_api_settings(self) -> None:
        """
        Validate all required settings are actually set.
        This is run on startup of the API which means that later on in the codebase
        we don't have to check for any non set values, we can just assume they are set.
        """
        required_fields = {
            "DIVBASE_ENV": self.api.environment,
            "FRONTEND_BASE_URL": self.api.frontend_base_url,
            "DATABASE_URL": self.database.url,
            "FLOWER_URL": self.flower.url,
            "FLOWER_USER": self.flower.user,
            "FLOWER_PASSWORD": self.flower.password,
            "JWT_SECRET_KEY": self.jwt.secret_key,
            "S3_ENDPOINT_URL": self.s3.endpoint_url,
            "S3_PRESIGNING_URL": self.s3.presigning_url,
            "S3_SERVICE_ACCOUNT_ACCESS_KEY": self.s3.access_key,
            "S3_SERVICE_ACCOUNT_SECRET_KEY": self.s3.secret_key,
        }
        for setting_name, setting in required_fields.items():
            if isinstance(setting, str) and setting == "NOT_SET":
                raise ValueError(f"A required environment variable was not set: {setting_name=}")
            if isinstance(setting, SecretStr):
                if setting.get_secret_value() == "NOT_SET":
                    raise ValueError(f"A required environment variable was not set: {setting_name=}")
                if self.api.environment not in ["local_dev", "test"] and setting.get_secret_value() == "badpassword":
                    raise ValueError(
                        f"A secret environment variable was set to badpassword for a non local environment: {setting_name=}"
                    )


# This instance can be imported and used throughout the application.
settings = Settings()
