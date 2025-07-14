"""
Functions that resolve for the CLI commands things like:
    - which bucket to use
    - which download directory to use
    - which DivBase API/S3 URL to use
Based on provider user input and their config file.
"""

from pathlib import Path

from divbase_tools.exceptions import BucketNameNotSpecifiedError
from divbase_tools.user_config import BucketConfig, load_user_config


def resolve_bucket(bucket_name: str | None, config_path: Path) -> BucketConfig:
    """
    Helper function to resolve the bucket to use for a CLI command.
    Falls back to the default bucket set in the user config if not explicitly provided.

    Once the bucket is resolved a BucketConfig object is returned,
    which contains name and URLs (S3+API) for the bucket.
    """
    config = load_user_config(config_path)
    if not bucket_name:
        bucket_name = config.default_bucket
    if not bucket_name:
        raise BucketNameNotSpecifiedError(config_path=config_path)
    return config.bucket_info(bucket_name)


def resolve_divbase_api_url(url: str | None, config_path: Path) -> str:
    """
    Helper function to resolve the DivBase API URL to use for a CLI command.

    If not provided by the user, it will take the default bucket's API URL from the user config.
    Otherwise raise an error.
    """
    if url:
        return url

    config = load_user_config(config_path=config_path)

    if not config.default_bucket:
        raise ValueError(
            "No default bucket is set in your user config. "
            "Please set a default bucket or specify the API URL using the --url option."
        )

    bucket_config = config.bucket_info(name=config.default_bucket)
    return bucket_config.divbase_url


def resolve_download_dir(download_dir: str | None, config_path: Path) -> Path:
    """
    Helper function to resolve the download directory to use for a CLI command involving downloading files.

    Priority given to `download_dir` argument, then if a default is set in the user config.
    Note: "." or None should default to the current working directory.
    """
    if not download_dir:
        config = load_user_config(config_path)
        download_dir = config.default_download_dir

    if download_dir and download_dir != ".":
        return Path(download_dir)
    return Path.cwd()
