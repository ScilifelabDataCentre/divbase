from pathlib import Path

import typer
from rich import print

from divbase_tools.cli_commands.config_resolver import resolve_project
from divbase_tools.cli_commands.user_config_cli import CONFIG_FILE_OPTION
from divbase_tools.services import (
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


@version_app.command("create")
def create_version(
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """Create a bucket versioning file that is stored inside the project's storage bucket."""
    project_config = resolve_project(project_name=project, config_path=config_file)
    create_version_object_command(project_config)
    print(f"Bucket versioning file created for project: '{project_config.name}'")


@version_app.command("add")
def add_version(
    name: str = typer.Argument(help="Name of the version (e.g., semantic version).", show_default=False),
    description: str = typer.Option(None, help="Optional description of the version."),
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """Add an entry to the bucket versioning file specfying the current state of all files in the project's storage bucket."""
    project_config = resolve_project(project_name=project, config_path=config_file)
    add_version_command(project_config=project_config, name=name, description=description)
    print(f"New version: '{name}' added to the project: '{project_config.name}'")


@version_app.command("list")
def list_versions(
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """List all entries in the bucket versioning file."""
    project_config = resolve_project(project_name=project, config_path=config_file)
    version_info = list_versions_command(project_config=project_config)

    if not version_info:
        print(f"No versions found for project: {project_config.name}.")
        return

    print("Project versions:")
    for version, details in version_info.items():
        desc = details["description"] or "No description provided"
        print(f"- '{version}': '{details['timestamp']}' : {desc}")


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
    deleted_version = delete_version_command(bucket_version=name, project_config=project_config)
    print(f"The version: '{deleted_version}' was deleted from the project: '{project_config.name}'")


@version_app.command("info")
def get_version_info(
    version: str = typer.Argument(help="Specific version to retrieve information for"),
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """Provide detailed information about a user specified project version, including all files present and their unique hashes."""
    project_config = resolve_project(project_name=project, config_path=config_file)
    files_at_version = list_files_at_version_command(project_config=project_config, bucket_version=version)

    if not files_at_version:
        print("No files were registered at this version.")
        return

    print(f"State of each file in the project: '{project_config.name}' at version: '{version}'")
    for object_name, hash in files_at_version.items():
        print(f"- '{object_name}' : '{hash}'")
