"""
Query subcommand for the DivBase CLI.

Submits queries (sample metadata and/or bcftools) to the DivBase API.

If sample metadata query:
    results are printed to the console.

If bcftools query:
    a task id is returned which can be used to check the status of the job.
    After task completed, a merged VCF file will be added to the project's storage bucket which can be downloaded by the user.


TODO:
- Ability to download results file given task id with the file cli?
-
"""

import logging
from pathlib import Path

import typer
from rich import print

from divbase_cli.cli_commands.shared_args_options import PROJECT_NAME_OPTION
from divbase_cli.cli_config import cli_settings
from divbase_cli.config_resolver import resolve_project
from divbase_cli.user_auth import make_authenticated_request
from divbase_lib.api_schemas.queries import (
    BcftoolsQueryRequest,
    SampleMetadataQueryRequest,
    SampleMetadataQueryTaskResult,
)

logger = logging.getLogger(__name__)


METADATA_TSV_ARGUMENT = typer.Option(
    cli_settings.METADATA_TSV_NAME, help="Name of the sample metadata TSV file in the project's data store on DivBase."
)

BCFTOOLS_ARGUMENT = typer.Option(
    ...,
    help="""
        String consisting of the bcftools command to run on the files returned by the tsv query.
        """,
)

# In 1 command this is required, other optional, hence only defining the text up here.
TSV_FILTER_HELP_TEXT = """String consisting of keys:values in the tsv file to filter on.
    The syntax is 'Key1:Value1,Value2;Key2:Value3,Value4', where the key
    are the column header names in the tsv, and values are the column values. 
    Multiple values for a key are separated by commas, and multiple keys are 
    separated by semicolons. When multple keys are provided, an intersect query 
    will be performed. E.g. 'Area:West of Ireland,Northern Portugal;Sex:F'.
    """


query_app = typer.Typer(
    help="Run queries on the VCF files stored in the project's data store on DivBase. Queries are run on the DivBase API",
    no_args_is_help=True,
)


@query_app.command("tsv")
def sample_metadata_query(
    filter: str = typer.Argument(
        ...,
        help=TSV_FILTER_HELP_TEXT,
    ),
    show_sample_results: bool = typer.Option(
        default=False,
        help="Print sample_ID and Filename results from the query.",
    ),
    metadata_tsv_name: str = METADATA_TSV_ARGUMENT,
    project: str | None = PROJECT_NAME_OPTION,
) -> None:
    """
    Query the tsv sidecar metadata file for the VCF files in the project's data store on DivBase.
    Returns the sample IDs and filenames that match the query.

    TODO: it perhaps be useful to set the default download_dir in the config so that we can
    look for files there? For now this code just uses file.parent as the download directory.
    TODO: handle when the name of the sample column is something other than Sample_ID
    """

    project_config = resolve_project(project_name=project)

    request_data = SampleMetadataQueryRequest(tsv_filter=filter, metadata_tsv_name=metadata_tsv_name)

    response = make_authenticated_request(
        method="POST",
        divbase_base_url=project_config.divbase_url,
        api_route=f"v1/query/sample-metadata/projects/{project_config.name}",
        json=request_data.model_dump(),
        timeout=20,  # This is longer than default (5), as api call response is query result, not a task-id.
    )

    results = SampleMetadataQueryTaskResult(**response.json())

    if results.warnings:
        print("[yellow]Warnings:[/yellow]")
        for warning in results.warnings:
            print(f"  • {warning}")
        print()

    if show_sample_results:
        print("[bright_blue]Name and file for each sample in query results:[/bright_blue]")
        for sample in results.sample_and_filename_subset:
            print(f"Sample ID: '{sample['Sample_ID']}', Filename: '{sample['Filename']}'")

    print(f"The results for the query ([bright_blue]{results.query_message}[/bright_blue]):")

    unique_sample_ids = results.unique_sample_ids or []
    unique_filenames = results.unique_filenames or []
    print(f"Unique Sample IDs: {unique_sample_ids}")
    print(f"Unique filenames: {unique_filenames}\n")


@query_app.command("bcftools-pipe")
def vcf_query(
    tsv_filter: str = typer.Option(None, help=TSV_FILTER_HELP_TEXT),
    samples: str | None = typer.Option(
        None,
        help="Comma-separated list of sample IDs. Mutually exclusive with --tsv-filter and --samples-file.",
    ),
    samples_file: Path | None = typer.Option(
        None,
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
        help=(
            "Path to a UTF-8 text file with one sample ID per line. Mutually exclusive with --tsv-filter and --samples."
        ),
    ),
    command: str = BCFTOOLS_ARGUMENT,
    metadata_tsv_name: str = METADATA_TSV_ARGUMENT,
    project: str | None = PROJECT_NAME_OPTION,
) -> None:
    """
    Submit a query to run on the DivBase API. A single, merged VCF file will be added to the project on success.

    TODO Error handling for subprocess calls.
    TODO: handle case empty results are returned from tsv_query()
    TODO what if the user just want to run bcftools on existing files in the bucket, without a tsv file query first?
    TODO what if a job fails and the user wants to re-run it? do we store temp files?
    TODO be consistent about input argument and options. when are they optional, how is that indicated in docstring? etc.
    TODO consider handling the bcftools command whitelist checks also on the CLI level since the error messages are nicer looking?
    TODO consider moving downloading of missing files elsewhere, since this is now done before the celery task
    """

    has_tsv_filter = tsv_filter is not None
    has_samples = samples is not None
    has_samples_file = samples_file is not None
    if sum([has_tsv_filter, has_samples, has_samples_file]) > 1:
        raise typer.BadParameter("Use only one of --tsv-filter, --samples, or --samples-file.")

    normalized_samples = _normalize_samples_input(samples=samples, samples_file=samples_file)

    project_config = resolve_project(project_name=project)

    request_data = BcftoolsQueryRequest(
        tsv_filter=tsv_filter,
        command=command,
        metadata_tsv_name=metadata_tsv_name,
        samples=normalized_samples,
    )

    response = make_authenticated_request(
        method="POST",
        divbase_base_url=project_config.divbase_url,
        api_route=f"v1/query/bcftools-pipe/projects/{project_config.name}",
        json=request_data.model_dump(),
    )

    task_id = response.json()
    print(f"Job submitted successfully with task id: {task_id}")


def _normalize_samples_input(samples: str | None, samples_file: Path | None) -> list[str] | None:
    """
    Normalize sample selection inputs from CLI options into a single list[str] or None.
    """
    if samples is not None:
        normalized = [sample.strip() for sample in samples.split(",") if sample.strip()]
        if not normalized:
            raise typer.BadParameter("--samples must contain at least one non-empty sample ID.")
        return normalized

    if samples_file is not None:
        normalized = [line.strip() for line in samples_file.read_text(encoding="utf-8").splitlines() if line.strip()]
        if not normalized:
            raise typer.BadParameter(f"Samples file is empty: {samples_file}")
        invalid_lines = [line for line in normalized if "," in line]
        if invalid_lines:
            raise typer.BadParameter(
                "Invalid --samples-file format: expected one sample ID per line with no commas. "
                "Use --samples for comma-separated input."
            )
        return normalized

    return None
