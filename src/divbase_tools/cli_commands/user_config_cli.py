"""
config subcommand for the divbase_tools CLI.

Controls the user config file stored at "~/.config/.divbase_tools.yaml" (unless specified otherwise).
"""

from pathlib import Path

import typer
from rich import print
from rich.console import Console
from rich.table import Table

from divbase_tools.user_config import (
    create_user_config,
    load_user_config,
)

DIVBASE_API_URL = "http://localhost:8000"
DIVBASE_S3_URL = "http://localhost:9000"
DEFAULT_CONFIG_PATH = Path.home() / ".config" / ".divbase_tools.yaml"

CONFIG_FILE_OPTION = typer.Option(
    DEFAULT_CONFIG_PATH,
    "--config",
    "-c",
    help="Path to your user configuration file. By default it is stored at ~/.config/.divbase_tools.yaml.",
)

config_app = typer.Typer(help="Manage your user configuration file for the DivBase CLI.", no_args_is_help=True)


@config_app.command("create")
def create_user_config_command(
    config_file: Path = typer.Option(DEFAULT_CONFIG_PATH, help="Where to store your config file locally on your pc."),
):
    """Create a user configuration file for the divbase_tools package."""
    create_user_config(config_path=config_file)
    print(f"User configuration file created at {config_file.resolve()}.")


@config_app.command("add-bucket")
def add_bucket_command(
    name: str = typer.Argument(..., help="Name of the storage bucket to add to the user configuration file."),
    divbase_url: str = typer.Option(DIVBASE_API_URL, help="DivBase API URL associated with this project."),
    s3_url: str = typer.Option(DIVBASE_S3_URL, help="S3 object store URL associated with this project."),
    make_default: bool = typer.Option(
        False,
        "--default",
        "-d",
        help="Set the bucket as the default bucket in the user configuration file.",
    ),
    config_file: Path = CONFIG_FILE_OPTION,
):
    """Add a storage bucket to the user configuration file."""
    config = load_user_config(config_file)
    bucket = config.add_bucket(
        name=name,
        divbase_url=divbase_url,
        s3_url=s3_url,
        is_default=make_default,
    )

    print(f"Bucket: '{bucket.name}' added to your config file located at {config_file.resolve()}.")
    print(f"DivBase URL: {bucket.divbase_url} and S3 URL: {bucket.s3_url} were set for this bucket.")

    if make_default:
        print(f"Bucket '{bucket.name}' is now set as your default bucket")
    else:
        print(f"To make '{bucket.name}' your default bucket you can run: 'divbase config set-default {bucket.name}'")


@config_app.command("remove-bucket")
def remove_bucket_command(
    name: str = typer.Argument(..., help="Name of the storage bucket to remove from the user configuration file."),
    config_file: Path = CONFIG_FILE_OPTION,
):
    """Remove a storage bucket from the user configuration file."""
    config = load_user_config(config_file)
    removed_bucket = config.remove_bucket(name)

    if not removed_bucket:
        print(f"The bucket '{name}' was not found in your config file located at {config_file.resolve()}.")
    else:
        print(f"The bucket '{removed_bucket}' was removed from your config.")


@config_app.command("set-default")
def set_default_bucket_command(
    name: str = typer.Argument(..., help="Name of the storage bucket to add to the user configuration file."),
    config_file: Path = CONFIG_FILE_OPTION,
):
    """Set the default bucket to use."""
    config = load_user_config(config_file)
    default_bucket_name = config.set_default_bucket(name=name)
    print(f"Default bucket is now set to '{default_bucket_name}'.")


@config_app.command("show-default")
def show_default_bucket_command(
    config_file: Path = CONFIG_FILE_OPTION,
) -> None:
    """Show the currently set default bucket."""
    config = load_user_config(config_file)

    if config.default_bucket:
        print(f"Default bucket is: '{config.default_bucket}'.")
    else:
        print("No default bucket is set in the user configuration file.")


@config_app.command("set-dload-dir")
def set_default_dload_dir_command(
    download_dir: str = typer.Argument(
        ...,
        help="""Set the default directory to download files to. 
        By default files are downloaded to the current working directory.
        You can specify an absolute path. 
        You can use '.' to refer to the directory you run the command from.""",
    ),
    config_file: Path = CONFIG_FILE_OPTION,
):
    """Set the default download dir"""
    config = load_user_config(config_file)
    dload_dir = config.set_default_download_dir(download_dir=download_dir)
    if dload_dir == ".":
        print("The default download directory will be whereever you run the command from.")
    else:
        print(f"The default download directory is now set to: {dload_dir}.")


@config_app.command("show")
def show_user_config(
    config_file: Path = CONFIG_FILE_OPTION,
):
    """Pretty print the contents of your current config file."""
    config = load_user_config(config_file)
    console = Console()

    table = Table(title=f"Your DivBase user configuration file's contents located at: '{config_file}'")
    table.add_column("Bucket Name", style="cyan")
    table.add_column("DivBase URL", style="green")
    table.add_column("S3 URL", style="green")
    table.add_column("Is default", style="yellow")

    if not config.default_download_dir:
        dload_dir_info = "Not specified, meaning the working directory of wherever you run the download command from."
    elif config.default_download_dir == ".":
        dload_dir_info = "Working directory of wherever you run the download command from."
    else:
        dload_dir_info = config.default_download_dir

    console.print(f"[bold]Default Download Directory:[/bold] {dload_dir_info}")

    if not config.buckets:
        console.print("[bold]No buckets defined in your user config file.[/bold]")
        console.print("You can add a bucket using the command: 'divbase config add-bucket <bucket_name>'")
        return

    for bucket in config.buckets:
        is_default = "Yes" if bucket.name == config.default_bucket else ""
        table.add_row(bucket.name, bucket.divbase_url, bucket.s3_url, is_default)
    console.print(table)
