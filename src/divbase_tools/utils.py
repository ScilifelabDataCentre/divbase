from pathlib import Path

from divbase_tools.exceptions import BucketNameNotSpecifiedError
from divbase_tools.user_config import load_user_config


def resolve_bucket_name(bucket_name: str | None, config_path: Path) -> str:
    """
    Helper function to resolve the bucket to use for a command.
    Falls back to the default bucket set in the user config if not explicitly provided.
    """
    if not bucket_name:
        config = load_user_config(config_path)
        bucket_name = config.default_bucket
    if not bucket_name:
        raise BucketNameNotSpecifiedError(config_path=config_path)

    return bucket_name
