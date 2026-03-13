import csv
import logging
from pathlib import Path

import typer
import yaml
from rich import print

from divbase_cli.cli_commands.shared_args_options import PROJECT_NAME_OPTION
from divbase_cli.config_resolver import resolve_project
from divbase_cli.user_auth import make_authenticated_request
from divbase_lib.api_schemas.vcf_dimensions import (
    DimensionsSamplesResult,
    DimensionsScaffoldsResult,
    DimensionsShowResult,
)
from divbase_lib.metadata_validator import SharedMetadataValidator

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
    sample_names_limit: int = typer.Option(
        20,
        "--sample-names-limit",
        min=1,
        help="Maximum number of sample names to display per list in terminal output.",
    ),
    sample_names_output: str | None = typer.Option(
        None,
        "--sample-names-output",
        help="Write full sample names to file instead of truncating in terminal output. "
        "Mutually exclusive with --sample-names-stdout.",
    ),
    sample_names_stdout: bool = typer.Option(
        False,
        "--sample-names-stdout",
        help="Print full sample names to stdout (useful for piping). Mutually exclusive with --sample-names-output.",
    ),
    project: str | None = PROJECT_NAME_OPTION,
) -> None:
    """
    Show the dimensions index file for a project.
    When running --unique-scaffolds, the sorting separates between numeric and non-numeric scaffold names.
    """

    project_config = resolve_project(project_name=project)

    # These two options are mutually exclusive. But due to how typer handles options, this error will only be raised if --sample-names-output has an input value (i.e. path).
    # If the path is missing, it will raise an error about the missing path argument instead...
    if sample_names_output and sample_names_stdout:
        typer.echo("Use only one of --sample-names-output or --sample-names-stdout.", err=True)
        raise typer.Exit(code=1)

    if unique_samples:
        response = make_authenticated_request(
            method="GET",
            divbase_base_url=project_config.divbase_url,
            api_route=f"v1/vcf-dimensions/projects/{project_config.name}/samples",
        )
        unique_sample_names_sorted = DimensionsSamplesResult(**response.json()).unique_samples
        sample_count = len(unique_sample_names_sorted)

        if sample_names_output:
            output_path = Path(sample_names_output)
            output_path.write_text("\n".join(unique_sample_names_sorted) + "\n")
            print(f"Wrote {sample_count} unique sample names to: {output_path}")
            return

        if sample_names_stdout:
            print("\n".join(unique_sample_names_sorted))
            return

        if sample_count > sample_names_limit:
            preview = unique_sample_names_sorted[:sample_names_limit]
            print(
                f"Unique sample names found across all the VCF files in the project (count: {sample_count}, showing first {sample_names_limit}):\n{preview}\n"
                f"To view all, use --sample-names-output <FILE> or --sample-names-stdout."
            )
            return

        print(
            f"Unique sample names found across all the VCF files in the project (count: {sample_count}):\n{unique_sample_names_sorted}"
        )
        return

    if unique_scaffolds:
        response = make_authenticated_request(
            method="GET",
            divbase_base_url=project_config.divbase_url,
            api_route=f"v1/vcf-dimensions/projects/{project_config.name}/scaffolds",
        )
        unique_scaffold_names_sorted = DimensionsScaffoldsResult(**response.json()).unique_scaffolds
        scaffold_count = len(unique_scaffold_names_sorted)
        print(
            f"Unique scaffold names found across all the VCF files in the project (count: {scaffold_count}):\n{unique_scaffold_names_sorted}"
        )
        return

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
            if sample_names_output or sample_names_stdout:
                _write_or_print_sample_names(
                    indexed_files=[record],
                    sample_names_output=sample_names_output,
                    sample_names_stdout=sample_names_stdout,
                )
                return

            _truncate_sample_names_in_entry(record, sample_names_limit)
            print(yaml.safe_dump(record, sort_keys=False))
        else:
            print(
                f"No entry found for filename: {filename}. Please check that the filename is correct and that it is a VCF file (extension: .vcf or .vcf.gz)."
                "\nHint: use 'divbase-cli files ls' to view all files in the project."
            )
        return

    if sample_names_output or sample_names_stdout:
        _write_or_print_sample_names(
            indexed_files=dimensions_info.get("indexed_files", []),
            sample_names_output=sample_names_output,
            sample_names_stdout=sample_names_stdout,
        )
        return

    for entry in dimensions_info.get("indexed_files", []):
        _truncate_sample_names_in_entry(entry, sample_names_limit)

    print(yaml.safe_dump(dimensions_info, sort_keys=False))


def _format_api_response_for_display_in_terminal(api_response: DimensionsShowResult) -> dict:
    """
    Convert the API response to a YAML-like format for display in the user's terminal.
    """

    def sort_scaffolds(scaffolds: list[str]) -> list[str]:
        numeric = sorted([int(s) for s in scaffolds if s.isdigit()])
        non_numeric = sorted([s for s in scaffolds if not s.isdigit()])
        return [str(n) for n in numeric] + non_numeric

    dimensions_list = []
    for entry in api_response.vcf_files:
        scaffolds = entry.get("scaffolds", [])
        sorted_scaffolds = sort_scaffolds(scaffolds)
        dimensions_entry = {
            "filename": entry["vcf_file_s3_key"],
            "file_version_ID_in_bucket": entry["s3_version_id"],
            "last_updated": entry.get("updated_at"),
            "dimensions": {
                "scaffold_count": len(sorted_scaffolds),
                "scaffolds": sorted_scaffolds,
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
    output_path: Path | None = typer.Option(
        None,
        "--output",
        "-o",
        file_okay=True,
        dir_okay=False,
        resolve_path=True,
        help="Path to the output TSV file to create. Defaults to sample_metadata_<project_name>.tsv in the current directory. If a file already exists at the given path, you will be prompted to confirm if you want to overwrite it.",
    ),
    project: str | None = PROJECT_NAME_OPTION,
) -> None:
    """
    Create a template sample metadata file (TSV format) pre-filled with the sample names from the project's VCF files based on the information stored in the project's VCF dimensions cache. Tip: run 'divbase-cli dimensions update' first to ensure that the VCF dimensions areup-to-date.
    """

    project_config = resolve_project(project_name=project)

    if output_path is None:
        output_path = Path.cwd() / f"sample_metadata_{project_config.name}.tsv"

    response = make_authenticated_request(
        method="GET",
        divbase_base_url=project_config.divbase_url,
        api_route=f"v1/vcf-dimensions/projects/{project_config.name}/samples",
    )
    unique_sample_names_sorted = DimensionsSamplesResult(**response.json()).unique_samples

    sample_count = len(unique_sample_names_sorted)
    print(
        f"There were {sample_count} unique samples found in the dimensions file for the {project_config.name} project."
    )

    if sample_count == 0:
        # Fallback in case there are no samples in the dimensions index. If no dimensions entry for the project
        # VCFDimensionsEntryMissingError will be returned. But for some reason, there are no samples in the VCF, this will catch that.
        print("No samples found for this project. No file written.")
        return

    # Check if file exists and prompt user for confirmation
    if output_path.exists():
        overwrite = typer.confirm(f"File '{output_path}' already exists. Do you want to overwrite it?")
        if not overwrite:
            print("File not written. Exiting.")
            return

    with open(output_path, mode="w", newline="") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        writer.writerow(["#Sample_ID"])
        for sample in unique_sample_names_sorted:
            writer.writerow([sample])

    print(f"A sample metadata template with these sample names was written to: {output_path}")

    # TODO perhaps add a message on how to fill in additional columns and how to upload the metadata file to DivBase?


@dimensions_app.command("validate-metadata-file")
def validate_metadata_template_versus_dimensions_and_formatting_constraints(
    input_path: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=False,
        resolve_path=True,
        help="Path to the input TSV file to validate.",
    ),
    untruncated: bool = typer.Option(
        False,
        "--untruncated",
        help="Show full (untruncated) validator lists (including sample mismatches and grouped warning row/value previews).",
    ),
    project: str | None = PROJECT_NAME_OPTION,
) -> None:
    """
    Validate a sample metadata TSV file (before you upload it to the project's data store) to check that it will work with DivBase queries.
    """
    # Client-side validation of a sidecar metadata TSV file, intended to be run before upload to DivBase.
    # Uses the SharedMetadataValidator (that is also used on the server-side) which checks for formatting errors and also validates that the sample names
    # in the TSV file match the sample names in the dimensions index for the project

    project_config = resolve_project(project_name=project)

    print(f"Validating local metadata file: {input_path}")
    print(f"Project: {project_config.name}\n")

    response = make_authenticated_request(
        method="GET",
        divbase_base_url=project_config.divbase_url,
        api_route=f"v1/vcf-dimensions/projects/{project_config.name}/samples",
    )
    unique_sample_names = DimensionsSamplesResult(**response.json()).unique_samples

    dimensions_sample_preview_limit = None if untruncated else 20

    shared_validator = SharedMetadataValidator(
        file_path=input_path,
        project_samples=set(unique_sample_names),
        skip_dimensions_check=False,
        dimensions_sample_preview_limit=dimensions_sample_preview_limit,
    )
    result = shared_validator.load_and_validate()

    errors = [error_entry.message for error_entry in result.errors]
    warnings = [warning_entry.message for warning_entry in result.warnings]
    stats = result.stats
    numeric_cols = result.numeric_columns
    string_cols = result.string_columns
    mixed_cols = result.mixed_type_columns

    print("[bold cyan]VALIDATION SUMMARY:[/bold cyan]")
    print(
        f"  Total columns: {getattr(stats, 'total_columns', 0)} ({getattr(stats, 'user_defined_columns', 0)} user-defined + 1 Sample_ID column)"
    )

    samples_in_tsv = getattr(stats, "samples_in_tsv", 0)
    samples_matching = getattr(stats, "samples_matching_project", 0)
    total_project = getattr(stats, "total_project_samples", 0)

    print(
        f"  Samples matching project VCF dimensions: {samples_matching}/{samples_in_tsv} (project has {total_project} total)"
    )

    print(f"  Numeric columns ({len(numeric_cols)}): {', '.join(numeric_cols) if numeric_cols else 'None'}")
    print(f"  String columns ({len(string_cols)}): {', '.join(string_cols) if string_cols else 'None'}")
    print(
        f"  Mixed-type columns treated as string ({len(mixed_cols)}): {', '.join(mixed_cols) if mixed_cols else 'None'}"
    )

    if getattr(stats, "has_multi_values", False):
        print("  Multi-value cells: Yes (Python list notation detected)")
    else:
        print("  Multi-value cells: No")

    empty_cells = getattr(stats, "empty_cells_per_column", {})
    if empty_cells:
        print(
            f"  User-defined columns with empty cells ({len(empty_cells)}): {', '.join(f'{col} ({count})' for col, count in empty_cells.items())}"
        )

    print()

    if errors:
        print("[red bold]ERRORS (must be fixed):[/red bold]")
        for error in errors:
            print(f"  - {error}")
        print()

    if warnings:
        print("[yellow bold]WARNINGS (should be reviewed):[/yellow bold]")
        for warning in warnings:
            print(f"  - {warning}")
        print()

    if not errors and not warnings:
        print(
            "[green bold]Validation passed![/green bold] The metadata file meets all DivBase requirements. The file is ready to be uploaded to your DivBase project with 'divbase-cli files upload'"
        )
    elif errors:
        print("[red bold]Validation failed![/red bold] Please fix the errors above before uploading.")
        raise typer.Exit(code=1)
    else:
        print("[yellow bold]Validation passed with warnings![/yellow bold] Review the warnings above.")

    # TODO: Add information about how to upload the validated metadata file to DivBase


def _truncate_sample_names_in_entry(entry: dict, sample_names_limit: int) -> None:
    """
    Truncate sample names output in terminal, while preserving full counts.
    """
    dimensions = entry.get("dimensions", {})
    sample_names = dimensions.get("sample_names", [])
    if len(sample_names) > sample_names_limit:
        dimensions["sample_names"] = sample_names[:sample_names_limit]
        dimensions["sample_names_note"] = (
            f"Showing first {sample_names_limit} of {len(sample_names)} samples. "
            "Use --sample-names-output <FILE> or --sample-names-stdout to view all."
        )


def _write_or_print_sample_names(
    indexed_files: list[dict],
    sample_names_output: str | None,
    sample_names_stdout: bool,
) -> None:
    """
    Export full sample names data from divbase-cli dimensions show
    either to a file or to stdout. If used without --unique-samples, it will also include the filename that each sample occurs in.
    """
    lines: list[str] = []
    for entry in indexed_files:
        file_name = entry.get("filename", "")
        sample_names = entry.get("dimensions", {}).get("sample_names", [])
        for sample_name in sample_names:
            lines.append(f"{file_name}\t{sample_name}")

    if sample_names_output:
        output_path = Path(sample_names_output)
        output_path.write_text("\n".join(lines) + ("\n" if lines else ""))
        print(f"Wrote {len(lines)} sample-name rows to: {output_path}")
        return

    if sample_names_stdout:
        print("\n".join(lines))
