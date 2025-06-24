"""
config subcommand for the divbase_tools CLI.

Controls the user config file stored at "~/.config/.divbase_tools.yaml" (unless specified otherwise).
"""

from pathlib import Path

import typer
from rich import print

from divbase_tools.user_config import (
    DEFAULT_CONFIG_PATH,
    add_bucket_to_config,
    create_user_config,
    get_default_bucket,
    load_user_config,
    remove_bucket_from_config,
    set_default_bucket,
)

config_app = typer.Typer(help="Manage your user configuration file for the DivBase CLI.", no_args_is_help=True)

# Needed by almost all commands, so we define it here.
CONFIG_FILE_OPTION = typer.Option(
    DEFAULT_CONFIG_PATH,
    "--config",
    "-c",
    help="Path to your user configuration file. By default it is stored at ~/.config/.divbase_tools.yaml.",
)


@config_app.command("create")
def create_user_config_command(
    config_file: Path = typer.Option(DEFAULT_CONFIG_PATH, help="Where to store your config file locally on your pc."),
):
    """Create a user configuration file for the divbase_tools package."""
    create_user_config(config_path=config_file)
    print(f"User configuration file created at {config_file.resolve()}.")


@config_app.command("add-bucket")
def add_bucket_to_config_command(
    bucket_name: str = typer.Argument(..., help="Name of the storage bucket to add to the user configuration file."),
    make_default: bool = typer.Option(
        False,
        "--default",
        "-d",
        help="Set the bucket as the default bucket in the user configuration file.",
    ),
    config_file: Path = CONFIG_FILE_OPTION,
):
    """Add a storage bucket to the user configuration file."""
    add_bucket_to_config(bucket_name=bucket_name, config_path=config_file, is_default=make_default)
    print(f"Bucket: '{bucket_name}' added to your config file located at {config_file.resolve()}.")


@config_app.command("remove-bucket")
def remove_bucket_command(
    bucket_name: str = typer.Argument(
        ..., help="Name of the storage bucket to remove from the user configuration file."
    ),
    config_file: Path = CONFIG_FILE_OPTION,
):
    """Remove a storage bucket from the user configuration file."""
    remove_bucket_from_config(bucket_name, config_path=config_file)
    print(f"Bucket '{bucket_name}' removed from your config file located at {config_file.resolve()}.")


@config_app.command("set-default")
def set_default_bucket_command(
    bucket_name: str = typer.Argument(..., help="Name of the storage bucket to add to the user configuration file."),
    config_file: Path = CONFIG_FILE_OPTION,
):
    """Set the default bucket to use."""
    set_default_bucket(bucket_name=bucket_name, config_path=config_file)
    print(f"Default bucket is now set to '{bucket_name}'.")


@config_app.command("show-default")
def show_default_bucket_command(
    config_file: Path = CONFIG_FILE_OPTION,
) -> None:
    """Show the currently set default bucket."""
    default = get_default_bucket(config_path=config_file)

    if default is None:
        print("No default bucket is set in the user configuration file.")
        return
    print(f"Default bucket is: '{default}'.")


@config_app.command("show")
def show_user_config(
    config_file: Path = CONFIG_FILE_OPTION,
):
    """Show the contents of your current config file."""
    config = load_user_config(config_file)
    print(f"Current user configuration from {config_file}:")
    for key, value in config.items():
        print(f"{key}: {value}")
