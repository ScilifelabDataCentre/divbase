from pathlib import Path

from divbase_tools.exceptions import BucketNameNotInConfigError, BucketNameNotSpecifiedError
from divbase_tools.user_config import BucketConfig, load_user_config


def resolve_bucket(bucket_name: str | None, config_path: Path) -> BucketConfig:
    """
    Helper function to resolve the bucket to use for a CLI command.
    Falls back to the default bucket set in the user config if not explicitly provided.

    Once the bucket is resovled a BucketConfig object is returned,
    which contains name and URLs (S3+API) for the bucket.
    """
    if not bucket_name:
        config = load_user_config(config_path)
        bucket_name = config.default_bucket
    if not bucket_name:
        raise BucketNameNotSpecifiedError(config_path=config_path)

    config = load_user_config(config_path)
    bucket = config.bucket_info(bucket_name)
    if not bucket:
        raise BucketNameNotInConfigError(config_path=config_path, bucket_name=bucket_name)
    return bucket


def resolve_download_dir(download_dir: str | None, config_path: Path) -> Path:
    """
    Helper function to resolve the download directory to use for a command that involves downloading files.

    Priority given to `download_dir` argument, then if a default is set in the user config.
    Note: "." or None should default to the current working directory.
    """
    if not download_dir:
        config = load_user_config(config_path)
        download_dir = config.default_download_dir

    if download_dir and download_dir != ".":
        return Path(download_dir)
    return Path.cwd()
