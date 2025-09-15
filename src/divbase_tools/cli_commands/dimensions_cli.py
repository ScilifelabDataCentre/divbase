import logging
from pathlib import Path

import httpx
import typer
import yaml

from divbase_tools.cli_commands.config_resolver import resolve_project
from divbase_tools.cli_commands.user_config_cli import CONFIG_FILE_OPTION
from divbase_tools.cli_commands.version_cli import PROJECT_NAME_OPTION
from divbase_tools.exceptions import VCFDimensionsFileMissingOrEmptyError
from divbase_tools.vcf_dimension_indexing import show_dimensions_command

logger = logging.getLogger(__name__)


dimensions_app = typer.Typer(
    no_args_is_help=True,
    help="Create and inspect dimensions (number of samples, number of variants, scaffold names) of the VCF files in a project",
)


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


@dimensions_app.command("show")
def show_dimensions_index(
    filename: str = typer.Option(
        None,
        "--filename",
        help="If set, will show only the entry for this VCF filename.",
    ),
    unique_scaffolds: bool = (
        typer.Option(
            False,
            "--unique-scaffolds",
            help="If set, will show all unique scaffold names found across all the VCF files in the project.",
        )
    ),
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
) -> None:
    """
    Show the dimensions index file for a project.
    When running --unique-scaffolds, the sorting separates between numeric and non-numeric scaffold names.
    """

    project_config = resolve_project(project_name=project, config_path=config_file)
    dimensions_info = show_dimensions_command(project_config=project_config)
    if not dimensions_info.get("dimensions"):
        raise VCFDimensionsFileMissingOrEmptyError(bucket_name=project_config.bucket_name)

    if filename:
        record = None
        for entry in dimensions_info.get("dimensions", []):
            if entry.get("filename") == filename:
                record = entry
                break
        if record:
            print(yaml.safe_dump(record, sort_keys=False))
        else:
            print(
                f"No entry found for filename: {filename}. Please check that the filename is correct and that it is a VCF file (extension: .vcf or .vcf.gz)."
                "\nHint: use 'divbase-cli files list' to view all files in the project."
            )
        return

    if unique_scaffolds:
        unique_scaffold_names = set()
        for entry in dimensions_info.get("dimensions", []):
            unique_scaffold_names.update(entry.get("dimensions", {}).get("scaffolds", []))

        numeric_scaffold_names = []
        non_numeric_scaffold_names = []
        for scaffold in unique_scaffold_names:
            if scaffold.isdigit():
                numeric_scaffold_names.append(int(scaffold))
            else:
                non_numeric_scaffold_names.append(scaffold)

        unique_scaffold_names_sorted = [str(n) for n in sorted(numeric_scaffold_names)] + sorted(
            non_numeric_scaffold_names
        )

        print(f"Unique scaffold names found across all the VCF files in the project:\n{unique_scaffold_names_sorted}")
        return

    print(yaml.safe_dump(dimensions_info, sort_keys=False))
