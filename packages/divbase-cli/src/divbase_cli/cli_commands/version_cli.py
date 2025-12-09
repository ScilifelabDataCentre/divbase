"""CLI commands for managing project versions in DivBase."""

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
    delete_version_command,
    get_version_details_command,
    list_versions_command,
)

PROJECT_NAME_OPTION = typer.Option(
    None,
    help="Name of the DivBase project, if not provided uses the default in your DivBase config file",
    show_default=False,
)

version_app = typer.Typer(
    no_args_is_help=True,
    help="Add, view and remove versions representing the state of all files in the entire project at the current timestamp.",
)


def format_timestamp(timestamp_str: str) -> str:
    """Format ISO timestamp to Europe/Stockholm format with timezone"""
    dt = datetime.fromisoformat(timestamp_str)
    cet_dt = dt.astimezone(ZoneInfo("Europe/Stockholm"))
    return cet_dt.strftime("%d/%m/%Y %H:%M:%S %Z")


@version_app.command("add")
def add_version(
    name: str = typer.Argument(help="Name of the version (e.g., semantic version).", show_default=False),
    description: str = typer.Option("", help="Optional description of the version."),
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """Add a new project version entry which specifies the current state of all files in the project at the current timestamp."""
    project_config = resolve_project(project_name=project, config_path=config_file)
    logged_in_url = ensure_logged_in(config_path=config_file, desired_url=project_config.divbase_url)

    add_version_response = add_version_command(
        name=name,
        description=description,
        project_name=project_config.name,
        divbase_base_url=logged_in_url,
    )
    print(f"New version: '{add_version_response.name}' added to the project: '{project_config.name}'")


@version_app.command("list")
def list_versions(
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """List all entries in the project versioning file."""
    project_config = resolve_project(project_name=project, config_path=config_file)
    logged_in_url = ensure_logged_in(config_path=config_file, desired_url=project_config.divbase_url)

    versions_info = list_versions_command(project_name=project_config.name, divbase_base_url=logged_in_url)

    if not versions_info:
        print(f"No versions found for project: {project_config.name}.")
        return

    console = Console()
    table = Table(title=f"Versions for {project_config.name}")
    table.add_column("Version", style="cyan", no_wrap=True)
    table.add_column("Created ", style="magenta")
    table.add_column("Description", style="green")

    for version in versions_info:
        name = version.name
        desc = version.description or "No description provided"
        created_at = format_timestamp(version.created_at)
        table.add_row(name, created_at, desc)

    console.print(table)


@version_app.command("info")
def get_version_info(
    version: str = typer.Argument(help="Specific version to retrieve information for"),
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """Provide detailed information about a user specified project version, including all files present and their unique hashes."""
    project_config = resolve_project(project_name=project, config_path=config_file)
    logged_in_url = ensure_logged_in(config_path=config_file, desired_url=project_config.divbase_url)

    files_at_version = get_version_details_command(
        project_name=project_config.name, divbase_base_url=logged_in_url, version_name=version
    )

    print(f"State of each file in the project: '{project_config.name}' at version: '{version}'")
    print(f"Created at: {format_timestamp(files_at_version.created_at)}")
    for object_name, hash in files_at_version.files.items():
        print(f"- '{object_name}' : '{hash}'")


@version_app.command("delete")
def delete_version(
    name: str = typer.Argument(help="Name of the version (e.g., semantic version).", show_default=False),
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """
    Delete a version entry in the project versioning table. This does not delete the files themselves.
    Deleted version entries older than 30 days will be permanently deleted.
    You can ask a DivBase admin to restore a deleted version within that time period.
    """
    project_config = resolve_project(project_name=project, config_path=config_file)
    logged_in_url = ensure_logged_in(config_path=config_file, desired_url=project_config.divbase_url)

    deleted_version = delete_version_command(
        project_name=project_config.name, divbase_base_url=logged_in_url, version_name=name
    )
    print(f"The version: '{deleted_version}' was deleted from the project: '{project_config.name}'")
