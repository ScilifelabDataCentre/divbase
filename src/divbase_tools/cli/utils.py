from pathlib import Path

import typer
from rich import print

from divbase_tools.user_config import get_default_bucket


def resolve_bucket_name(bucket_name: str | None, config_path: Path) -> str:
    """
    Helper function to resolve the bucket to use for a command.
    Fallsback to the default bucket in user config if not explicitly provided.
    """
    if not bucket_name:
        bucket_name = get_default_bucket(config_path)
    if not bucket_name:
        print(
            "No bucket name provided. Please provide a bucket name or set a default bucket in your user configuration file."
        )
        raise typer.Exit(code=1)

    return bucket_name
