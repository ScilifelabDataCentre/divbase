import logging
from pathlib import Path

import httpx
import typer
import yaml

from divbase_cli.cli_commands.user_config_cli import CONFIG_FILE_OPTION
from divbase_cli.cli_commands.version_cli import PROJECT_NAME_OPTION
from divbase_cli.config_resolver import resolve_project
from divbase_cli.user_auth import make_authenticated_request
from divbase_lib.exceptions import VCFDimensionsFileMissingOrEmptyError

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
    response = httpx.post(f"{project_config.divbase_url}/v1/dimensions/update/", params=params)
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

    try:
        response = make_authenticated_request(
            method="GET",
            divbase_base_url=project_config.divbase_url,
            api_route=f"v1/vcf-dimensions/list/user-project-name/{project_config.name}",
        )
        vcf_dimensions_data = response.json()
    except httpx.HTTPStatusError as e:
        if e.response.status_code == 404:
            raise VCFDimensionsFileMissingOrEmptyError(bucket_name=project_config.name) from e
        elif e.response.status_code == 403:
            print(f"Error: You don't have access to project '{project_config.name}'")
            raise typer.Exit(1) from None
        raise

    if not vcf_dimensions_data.get("vcf_files"):
        raise VCFDimensionsFileMissingOrEmptyError(bucket_name=project_config.name)

    if not vcf_dimensions_data.get("vcf_files"):
        raise VCFDimensionsFileMissingOrEmptyError(bucket_name=project_config.name)

    dimensions_info = _format_api_response_for_display_in_terminal(vcf_dimensions_data)

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


def _format_api_response_for_display_in_terminal(api_response: dict) -> dict:
    """
    Convert the new API response format to a YAML-like format.
    """
    dimensions_list = []

    for entry in api_response.get("vcf_files", []):
        dimensions_entry = {
            "filename": entry["vcf_file_s3_key"],
            "file_version_ID_in_bucket": entry["s3_version_id"],
            "last_updated": entry.get("indexed_at") or entry.get("updated_at"),
            "dimensions": {
                "scaffolds": entry.get("scaffolds", []),
                "sample_count": entry.get("sample_count", 0),
                "sample_names": entry.get("samples", []),
                "variants": entry.get("variant_count", 0),
            },
        }
        dimensions_list.append(dimensions_entry)

    skipped_list = []
    for entry in api_response.get("skipped_files", []):
        skipped_entry = {
            "filename": entry["vcf_file_s3_key"],
            "file_version_ID_in_bucket": entry["s3_version_id"],
            "skip_reason": entry.get("skip_reason", "unknown"),
        }
        skipped_list.append(skipped_entry)

    return {
        "indexed_files": dimensions_list,
        "skipped_files": skipped_list,
    }
