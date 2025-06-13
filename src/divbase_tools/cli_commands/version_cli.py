from pathlib import Path

import typer
from rich import print

from divbase_tools.cli_commands.user_config_cli import CONFIG_PATH_OPTION
from divbase_tools.services import (
    add_version_command,
    create_version_object_command,
    delete_version_command,
    list_versions_command,
)
from divbase_tools.utils import resolve_bucket_name

BUCKET_NAME_OPTION = typer.Option(None, help="Name of the storage bucket for the project.", show_default=False)

version_app = typer.Typer(
    no_args_is_help=True, help="Version the state of all files in the entire bucket at a given timestamp."
)


@version_app.command("create")
def create_version(
    bucket_name: str | None = BUCKET_NAME_OPTION,
    config_path: Path = CONFIG_PATH_OPTION,
):
    """Create a bucket versioning file that is stored in the bucket."""
    bucket_name = resolve_bucket_name(bucket_name=bucket_name, config_path=config_path)
    create_version_object_command(bucket_name=bucket_name, config_path=config_path)
    print(f"Bucket versioning file created in bucket: '{bucket_name}'")


@version_app.command("add")
def add_version(
    name: str = typer.Argument(help="Name of the version (e.g., semantic version).", show_default=False),
    description: str = typer.Option(None, help="Optional description of the version."),
    bucket_name: str | None = BUCKET_NAME_OPTION,
    config_path: Path = CONFIG_PATH_OPTION,
):
    """Add an entry to the bucket versioning file specfying the current state of all files in the bucket."""
    bucket_name = resolve_bucket_name(bucket_name=bucket_name, config_path=config_path)
    add_version_command(bucket_name, name, description, config_path=config_path)
    print(f"New version: '{name}' added to the bucket: '{bucket_name}'")


@version_app.command("list")
def list_versions(
    bucket_name: str | None = BUCKET_NAME_OPTION,
    config_path: Path = CONFIG_PATH_OPTION,
):
    """List all entries in the bucket versioning file."""
    bucket_name = resolve_bucket_name(bucket_name=bucket_name, config_path=config_path)
    version_info = list_versions_command(bucket_name, config_path=config_path)

    if not version_info:
        print(f"No versions found for bucket: {bucket_name}.")
        return

    print("Bucket versions:")
    for version, details in version_info.items():
        desc = details["description"] or "No description provided"
        print(f"- '{version}': '{details['timestamp']}' : {desc}")


@version_app.command("delete")
def delete_version(
    name: str = typer.Argument(help="Name of the version (e.g., semantic version).", show_default=False),
    bucket_name: str | None = BUCKET_NAME_OPTION,
    config_path: Path = CONFIG_PATH_OPTION,
):
    """Delete an entry in the bucket versioning file specfying a specific state of all files in the bucket. Does not delete the files themselves."""
    bucket_name = resolve_bucket_name(bucket_name=bucket_name, config_path=config_path)
    deleted_version = delete_version_command(bucket_version=name, bucket_name=bucket_name, config_path=config_path)
    print(f"The version: '{deleted_version}' was deleted from the bucket: '{bucket_name}'")
