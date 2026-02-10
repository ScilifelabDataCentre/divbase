import csv
import logging
from pathlib import Path

import typer
import yaml

from divbase_cli.cli_commands.shared_args_options import PROJECT_NAME_OPTION
from divbase_cli.config_resolver import resolve_project
from divbase_cli.user_auth import make_authenticated_request
from divbase_lib.api_schemas.vcf_dimensions import DimensionsShowResult

logger = logging.getLogger(__name__)


dimensions_app = typer.Typer(
    no_args_is_help=True,
    help="Create and inspect dimensions (number of samples, number of variants, scaffold names) of the VCF files in a project",
)


@dimensions_app.command("update")
def update_dimensions_index(
    project: str | None = PROJECT_NAME_OPTION,
) -> None:
    """Calculate and add the dimensions of a VCF file to the dimensions index file in the project."""

    project_config = resolve_project(project_name=project)

    response = make_authenticated_request(
        method="PUT",
        divbase_base_url=project_config.divbase_url,
        api_route=f"v1/vcf-dimensions/projects/{project_config.name}",
    )

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
    unique_samples: bool = (
        typer.Option(
            False,
            "--unique-samples",
            help="If set, will show all unique sample names found across all the VCF files in the project.",
        )
    ),
    project: str | None = PROJECT_NAME_OPTION,
) -> None:
    """
    Show the dimensions index file for a project.
    When running --unique-scaffolds, the sorting separates between numeric and non-numeric scaffold names.
    """

    project_config = resolve_project(project_name=project)

    response = make_authenticated_request(
        method="GET",
        divbase_base_url=project_config.divbase_url,
        api_route=f"v1/vcf-dimensions/projects/{project_config.name}",
    )
    vcf_dimensions_data = DimensionsShowResult(**response.json())

    dimensions_info = _format_api_response_for_display_in_terminal(vcf_dimensions_data)

    if filename:
        record = None
        for entry in dimensions_info.get("indexed_files", []):
            if entry.get("filename") == filename:
                record = entry
                break
        if record:
            print(yaml.safe_dump(record, sort_keys=False))
        else:
            print(
                f"No entry found for filename: {filename}. Please check that the filename is correct and that it is a VCF file (extension: .vcf or .vcf.gz)."
                "\nHint: use 'divbase-cli files ls' to view all files in the project."
            )
        return

    if unique_scaffolds:
        # TODO for scalability: implement this as a separate CRUD instead of parsing all data on the client side
        unique_scaffold_names = set()
        for entry in dimensions_info.get("indexed_files", []):
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

    if unique_samples:
        # TODO for scalability: implement this as a separate CRUD instead of parsing all data on the client side
        unique_sample_names = set()
        for entry in dimensions_info.get("indexed_files", []):
            unique_sample_names.update(entry.get("dimensions", {}).get("sample_names", []))

        unique_sample_names_sorted = sorted(unique_sample_names)

        print(f"Unique sample names found across all the VCF files in the project:\n{unique_sample_names_sorted}")
        return

    print(yaml.safe_dump(dimensions_info, sort_keys=False))


def _format_api_response_for_display_in_terminal(api_response: DimensionsShowResult) -> dict:
    """
    Convert the API response to a YAML-like format for display in the user's terminal.
    """
    dimensions_list = []
    for entry in api_response.vcf_files:
        dimensions_entry = {
            "filename": entry["vcf_file_s3_key"],
            "file_version_ID_in_bucket": entry["s3_version_id"],
            "last_updated": entry.get("updated_at"),
            "dimensions": {
                "scaffolds": entry.get("scaffolds", []),
                "sample_count": entry.get("sample_count", 0),
                "sample_names": entry.get("samples", []),
                "variants": entry.get("variant_count", 0),
            },
        }
        dimensions_list.append(dimensions_entry)

    skipped_list = []
    for entry in api_response.skipped_files:
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


@dimensions_app.command("create-metadata-template")
def create_metadata_template_with_project_samples_names(
    project: str | None = PROJECT_NAME_OPTION,
) -> None:
    """
    Use the samples index in a projects dimensions cache to create a TSV metadata template file
    that has the sample names as pre-filled as the first column.
    """

    # TODO this duplicates some code with show_dimensions_index() above. A refactoring should probably include creating a separate CRUD function
    # so that the client does not need to parse all data.

    project_config = resolve_project(project_name=project)

    response = make_authenticated_request(
        method="GET",
        divbase_base_url=project_config.divbase_url,
        api_route=f"v1/vcf-dimensions/projects/{project_config.name}",
    )
    vcf_dimensions_data = DimensionsShowResult(**response.json())

    dimensions_info = _format_api_response_for_display_in_terminal(vcf_dimensions_data)

    unique_sample_names = set()
    for entry in dimensions_info.get("indexed_files", []):
        unique_sample_names.update(entry.get("dimensions", {}).get("sample_names", []))

    unique_sample_names_sorted = sorted(unique_sample_names)
    sample_count = len(unique_sample_names_sorted)
    print(
        f"There were {sample_count} unique samples found in the dimensions file for the {project_config.name} project."
    )

    if sample_count == 0:
        # Fallback in case there are no samples in the dimensions index. If no dimensions entry for the project
        # VCFDimensionsEntryMissingError will be returned. But for some reason, there are no samples in the VCF, this will catch that.
        print("No samples found for this project. No file written.")
        return

    output_filename = "divbase_metadata_template.tsv"
    output_path = Path.cwd() / output_filename

    with open(output_path, mode="w", newline="") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        writer.writerow(["#Sample_ID"])
        for sample in unique_sample_names_sorted:
            writer.writerow([sample])

    print(f"A sample metadata template with these sample names was written to: {output_path}")

    # TODO perhaps add a message on how to fill in additional columns and how to upload the metadata file to DivBase?
