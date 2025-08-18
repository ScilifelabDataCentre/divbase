from pathlib import Path

import httpx
import typer

from divbase_tools.cli_commands.config_resolver import resolve_project
from divbase_tools.cli_commands.user_config_cli import CONFIG_FILE_OPTION
from divbase_tools.s3_client import create_s3_file_manager
from divbase_tools.vcf_dimension_indexing import VCFDimensionIndexManager

PROJECT_NAME_OPTION = typer.Option(
    None,
    help="Name of the DivBase project, if not provided uses the default in your DivBase config file",
    show_default=False,
)


dimensions_app = typer.Typer(
    no_args_is_help=True,
    help="Create and inspect dimensions (number of samples, number of variants, scaffold names) of the VCF files in a project",
)


@dimensions_app.command("create")
def create_dimensions_index(
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """Create a dimensions indexing file that is stored inside the project's storage bucket."""
    project_config = resolve_project(project_name=project, config_path=config_file)
    s3_file_manager = create_s3_file_manager(project_config.s3_url)
    VCFDimensionIndexManager(bucket_name=project_config.bucket_name, s3_file_manager=s3_file_manager)
    print(f"Dimensions indexing file created for project: '{project_config.name}'")

    # TODO this should check if there is a file already, and if there is it should do nothing and ask the user to run other commands instead
    # TODO submit a celery task to calculate all dimensions for VCF files


@dimensions_app.command("update")
def update_dimensions_index(
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
) -> None:
    """Calculate and add the dimensions of a VCF file to the dimensions index file in the project."""

    project_config = resolve_project(project_name=project, config_path=config_file)

    params = {
        "project": project_config.name,
    }
    response = httpx.post(f"{project_config.divbase_url}/dimensions/update/", params=params)
    response.raise_for_status()

    task_id = response.json()
    print(f"Job submitted successfully with task id: {task_id}")


# TODO add command to list all scaffolds available in bucket
