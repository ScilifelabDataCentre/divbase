from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import typer
from rich.console import Console
from rich.table import Table

from divbase_cli.cli_commands.user_config_cli import CONFIG_FILE_OPTION
from divbase_cli.config_resolver import ensure_logged_in, resolve_project
from divbase_cli.services import (
    add_version_command,
    create_version_object_command,
    delete_version_command,
    list_files_at_version_command,
    list_versions_command,
)

PROJECT_NAME_OPTION = typer.Option(
    None,
    help="Name of the DivBase project, if not provided uses the default in your DivBase config file",
    show_default=False,
)

version_app = typer.Typer(
    no_args_is_help=True,
    help="Version the state of all files in the entire projects storage bucket at a given timestamp.",
)


def format_timestamp(timestamp_str: str) -> str:
    """Format ISO timestamp to Europe/Stockholm format with timezone"""
    dt = datetime.fromisoformat(timestamp_str)
    cet_dt = dt.astimezone(ZoneInfo("Europe/Stockholm"))
    return cet_dt.strftime("%d/%m/%Y %H:%M:%S %Z")


@version_app.command("create")
def create_version(
    name: str = typer.Option(default="v0.0.0", help="Name of the version (e.g., semantic version)."),
    description: str = typer.Option("", help="Optional description of the version."),
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """Create a bucket versioning file that is stored inside the project's storage bucket."""
    project_config = resolve_project(project_name=project, config_path=config_file)
    logged_in_url = ensure_logged_in(config_path=config_file, desired_url=project_config.divbase_url)

    new_version = create_version_object_command(
        project_name=project_config.name,
        divbase_base_url=logged_in_url,
        version_name=name,
        description=description if description else "",
    )
    print(
        f"Bucket versioning file created for project: '{project_config.name}'\n"
        f"with initial version named: '{new_version.name}'\n"
        f" and description: '{new_version.description}'\n"
    )


@version_app.command("add")
def add_version(
    name: str = typer.Argument(help="Name of the version (e.g., semantic version).", show_default=False),
    description: str = typer.Option("", help="Optional description of the version."),
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """Add an entry to the bucket versioning file specfying the current state of all files in the project's storage bucket."""
    project_config = resolve_project(project_name=project, config_path=config_file)
    logged_in_url = ensure_logged_in(config_path=config_file, desired_url=project_config.divbase_url)

    add_version_command(
        name=name,
        description=description,
        project_name=project_config.name,
        divbase_base_url=logged_in_url,
    )
    print(f"New version: '{name}' added to the project: '{project_config.name}'")


@version_app.command("list")
def list_versions(
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """List all entries in the bucket versioning file."""
    project_config = resolve_project(project_name=project, config_path=config_file)
    logged_in_url = ensure_logged_in(config_path=config_file, desired_url=project_config.divbase_url)

    version_info = list_versions_command(project_name=project_config.name, divbase_base_url=logged_in_url)

    if not version_info:
        print(f"No versions found for project: {project_config.name}.")
        return

    console = Console()
    table = Table(title=f"Versions for {project_config.name}")
    table.add_column("Version", style="cyan", no_wrap=True)
    table.add_column("Created ", style="magenta")
    table.add_column("Description", style="green")

    for version, details in version_info.items():
        desc = details.description or "No description provided"
        formatted_timestamp = format_timestamp(details.timestamp)
        table.add_row(version, formatted_timestamp, desc)

    console.print(table)


@version_app.command("delete")
def delete_version(
    name: str = typer.Argument(help="Name of the version (e.g., semantic version).", show_default=False),
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """
    Delete an entry in the bucket versioning file specfying a specific state of all files in the project's storage bucket.
    Does not delete the files themselves.
    """
    project_config = resolve_project(project_name=project, config_path=config_file)
    logged_in_url = ensure_logged_in(config_path=config_file, desired_url=project_config.divbase_url)

    deleted_version = delete_version_command(
        project_name=project_config.name, divbase_base_url=logged_in_url, version_name=name
    )
    print(f"The version: '{deleted_version}' was deleted from the project: '{project_config.name}'")


@version_app.command("info")
def get_version_info(
    version: str = typer.Argument(help="Specific version to retrieve information for"),
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """Provide detailed information about a user specified project version, including all files present and their unique hashes."""
    project_config = resolve_project(project_name=project, config_path=config_file)
    logged_in_url = ensure_logged_in(config_path=config_file, desired_url=project_config.divbase_url)

    files_at_version = list_files_at_version_command(
        project_name=project_config.name, divbase_base_url=logged_in_url, bucket_version=version
    )

    if not files_at_version:
        print("No files were registered at this version.")
        return

    print(f"State of each file in the project: '{project_config.name}' at version: '{version}'")
    for object_name, hash in files_at_version.items():
        print(f"- '{object_name}' : '{hash}'")
