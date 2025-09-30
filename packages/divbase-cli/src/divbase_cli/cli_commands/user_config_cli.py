"""
config subcommand for the divbase-cli package.

Controls the user config file stored at "~/.config/.divbase_tools.yaml" (unless specified otherwise).
"""

from pathlib import Path

import typer
from rich import print
from rich.console import Console
from rich.table import Table

from divbase_cli.config import settings
from divbase_cli.user_config import (
    create_user_config,
    load_user_config,
)

CONFIG_FILE_OPTION = typer.Option(
    settings.DEFAULT_CONFIG_PATH,
    "--config",
    "-c",
    help="Path to your user configuration file. By default it is stored at ~/.config/.divbase_tools.yaml.",
)

config_app = typer.Typer(help="Manage your user configuration file for the DivBase CLI.", no_args_is_help=True)


@config_app.command("create")
def create_user_config_command(
    config_file: Path = typer.Option(
        settings.DEFAULT_CONFIG_PATH, help="Where to store your config file locally on your pc."
    ),
):
    """Create a user configuration file for the divbase-cli tool."""
    create_user_config(config_path=config_file)
    print(f"User configuration file created at {config_file.resolve()}.")


@config_app.command("add-project")
def add_project_command(
    name: str = typer.Argument(..., help="Name of the project to add to the user configuration file."),
    divbase_url: str = typer.Option(
        settings.DEFAULT_DIVBASE_API_URL, help="DivBase API URL associated with this project."
    ),
    s3_url: str = typer.Option(settings.DEFAULT_S3_URL, help="S3 object store URL associated with this project."),
    make_default: bool = typer.Option(
        False,
        "--default",
        "-d",
        help="Set this project as the default project in the user configuration file.",
    ),
    config_file: Path = CONFIG_FILE_OPTION,
):
    """Add a new project to the user configuration file."""
    config = load_user_config(config_file)
    project = config.add_project(
        name=name,
        divbase_url=divbase_url,
        s3_url=s3_url,
        is_default=make_default,
    )

    print(f"Project: '{project.name}' added to your config file located at {config_file.resolve()}.")
    print(f"DivBase URL: {project.divbase_url} and S3 URL: {project.s3_url} were set for this project.")

    if make_default:
        print(f"Project '{project.name}' is now set as your default project")
    else:
        print(f"To make '{project.name}' your default project you can run: 'divbase config set-default {project.name}'")


@config_app.command("remove-project")
def remove_project_command(
    name: str = typer.Argument(..., help="Name of the project to remove from the user configuration file."),
    config_file: Path = CONFIG_FILE_OPTION,
):
    """Remove a project from the user configuration file."""
    config = load_user_config(config_file)
    removed_project = config.remove_project(name)

    if not removed_project:
        print(f"The project '{name}' was not found in your config file located at {config_file.resolve()}.")
    else:
        print(f"The project '{removed_project}' was removed from your config.")


@config_app.command("set-default")
def set_default_project_command(
    name: str = typer.Argument(..., help="Name of the project to add to the user configuration file."),
    config_file: Path = CONFIG_FILE_OPTION,
):
    """Set the default project to use."""
    config = load_user_config(config_file)
    default_project_name = config.set_default_project(name=name)
    print(f"Default project is now set to '{default_project_name}'.")


@config_app.command("show-default")
def show_default_project_command(
    config_file: Path = CONFIG_FILE_OPTION,
) -> None:
    """Show the currently set default project."""
    config = load_user_config(config_file)

    if config.default_project:
        print(f"Default project is: '{config.default_project}'.")
    else:
        print("No default project is set in the user configuration file.")


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

    if not config.projects:
        console.print("[bold]No projects defined in your user config file.[/bold]")
        console.print("You can add a project using the command: 'divbase config add-project <project_name>'")
        return

    for project in config.projects:
        is_default = "Yes" if project.name == config.default_project else ""
        table.add_row(project.name, project.divbase_url, project.s3_url, is_default)
    console.print(table)
