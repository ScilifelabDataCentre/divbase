"""
Config for the DivBase celery workers

This module creates a single instance of the WorkerSettings class,
which other modules can import directly to access configuration settings.

This settings object is validated in the worker process init function.

For local development these settings are defined by environment variables,
which are set in the docker compose file.
"""

import os
from dataclasses import dataclass, field

from pydantic import SecretStr


@dataclass
class WorkerGeneralSettings:
    """General configuration settings for the worker."""

    environment: str = os.getenv("DIVBASE_ENV", "NOT_SET")
    sync_url: SecretStr = SecretStr(os.getenv("SYNC_DATABASE_URL", "NOT_SET"))
    broker_url: SecretStr = SecretStr(os.getenv("CELERY_BROKER_URL", "NOT_SET"))
    result_backend: SecretStr = SecretStr(os.getenv("CELERY_RESULT_BACKEND", "NOT_SET"))


@dataclass
class WorkerS3Settings:
    """S3 configuration settings for the worker."""

    endpoint_url: str = os.getenv("S3_ENDPOINT_URL", "http://host.docker.internal:9000")
    bucket_prefix: str = os.getenv("S3_BUCKET_PREFIX", "NOT_SET")
    access_key: SecretStr = SecretStr(os.getenv("S3_SERVICE_ACCOUNT_ACCESS_KEY", "NOT_SET"))
    secret_key: SecretStr = SecretStr(os.getenv("S3_SERVICE_ACCOUNT_SECRET_KEY", "NOT_SET"))


@dataclass
class WorkerMetricsSettings:
    """Worker metrics configuration settings."""

    # ENABLE_WORKER_METRICS controls whether the Prometheus metrics server is started (system metrics, etc.)
    enabled: bool = os.getenv("ENABLE_WORKER_METRICS", "1") == "1"
    # ENABLE_WORKER_METRICS_PER_TASK controls whether per-task metrics (task/bcftools/VCF download) are collected and exposed
    enabled_per_task: bool = os.getenv("ENABLE_WORKER_METRICS_PER_TASK", "1") == "1"
    # Prometheus scrapes every 15 seconds in DivBase setup. A TTL of 5 min means it is available for 20 scrapes.
    # Once Prometheus has scraped it, it will store the data in its own volume for its retention time (default 15d).
    cache_ttl_minutes: int = int(os.getenv("TASK_METRICS_CACHE_TTL_MINUTES", "5"))


@dataclass
class WorkerCronSettings:
    """
    Cron task retention settings for the worker.
    Used to define max age of items in the db or S3 before automatic hard/soft deletion by cron tasks.
    """

    stuck_pending_hours: int = int(os.getenv("STUCK_PENDING_STATUS_HOURS", "168"))  # 168 h = 7 days
    stuck_started_hours: int = int(os.getenv("STUCK_STARTED_STATUS_HOURS", "168"))  # 168 h = 7 days
    task_retention_days: int = int(os.getenv("TASK_RETENTION_DAYS", "30"))
    revoked_token_retention_days: int = int(os.getenv("REVOKED_TOKEN_RETENTION_DAYS", "7"))
    soft_deleted_files_retention_days: int = int(os.getenv("SOFT_DELETED_FILES_RETENTION_DAYS", "30"))
    soft_deleted_project_version_retention_days: int = int(
        os.getenv("SOFT_DELETED_PROJECT_VERSION_RETENTION_DAYS", "30")
    )


@dataclass
class WorkerSettings:
    """Configuration settings for DivBase Celery workers."""

    general: WorkerGeneralSettings = field(default_factory=WorkerGeneralSettings)
    s3: WorkerS3Settings = field(default_factory=WorkerS3Settings)
    metrics: WorkerMetricsSettings = field(default_factory=WorkerMetricsSettings)
    cron: WorkerCronSettings = field(default_factory=WorkerCronSettings)

    def validate(self) -> None:
        """Validate all required settings are set. Called on worker process startup."""
        required_fields = {
            "DIVBASE_ENV": self.general.environment,
            "SYNC_DATABASE_URL": self.general.sync_url,
            "CELERY_BROKER_URL": self.general.broker_url,
            "CELERY_RESULT_BACKEND": self.general.result_backend,
            "S3_SERVICE_ACCOUNT_ACCESS_KEY": self.s3.access_key,
            "S3_SERVICE_ACCOUNT_SECRET_KEY": self.s3.secret_key,
            "S3_BUCKET_PREFIX": self.s3.bucket_prefix,
        }
        for setting_name, setting in required_fields.items():
            if isinstance(setting, str) and setting == "NOT_SET":
                raise ValueError(f"A required environment variable was not set: {setting_name=}")
            if isinstance(setting, SecretStr):
                if setting.get_secret_value() == "NOT_SET":
                    raise ValueError(f"A required environment variable was not set: {setting_name=}")
                if (
                    self.general.environment not in ["local_dev", "test"]
                    and setting.get_secret_value() == "badpassword"
                ):
                    raise ValueError(
                        f"A secret environment variable was set to badpassword for a non local environment: {setting_name=}"
                    )

        if self.metrics.enabled_per_task and not self.metrics.enabled:
            raise ValueError(
                "ENABLE_WORKER_METRICS_PER_TASK cannot be set if ENABLE_WORKER_METRICS is not set."
                "Set both to '1' to enable per-task metrics collection, or set ENABLE_WORKER_METRICS_PER_TASK to '0' to disable per-task metrics collection."
            )


# This instance can be imported and used across the worker codebase to access settings.
worker_settings = WorkerSettings()
