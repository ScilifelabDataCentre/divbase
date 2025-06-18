import datetime
from pathlib import Path

from divbase_tools.exceptions import BucketNameNotSpecifiedError
from divbase_tools.user_config import get_default_bucket


def resolve_bucket_name(bucket_name: str | None, config_path: Path) -> str:
    """
    Helper function to resolve the bucket to use for a command.
    Falls back to the default bucket set in the user config if not explicitly provided.
    """
    if not bucket_name:
        bucket_name = get_default_bucket(config_path)
    if not bucket_name:
        raise BucketNameNotSpecifiedError(config_path=config_path)

    return bucket_name


def format_unix_timestamp(timestamp):
    """The flower task status API returns timestamps as integers or floats.
    This function formats them into a human-readable string"""

    # TODO check how this handles timezones
    try:
        if isinstance(timestamp, (int, float)):
            dt = datetime.datetime.fromtimestamp(timestamp)
            return dt.strftime("%Y-%m-%d %H:%M:%S")
        return str(timestamp)
    except Exception:
        return str(timestamp)
