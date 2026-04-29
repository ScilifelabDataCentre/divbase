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

import stamina
import typer
from rich import print

from divbase_cli.cli_commands.shared_args_options import PROJECT_NAME_OPTION
from divbase_cli.cli_config import cli_settings
from divbase_cli.config_resolver import ensure_logged_in, resolve_project, resolve_url_for_non_project_specific_commands
from divbase_cli.retries import (
    retry_polling_until_final_or_retryable_api_errors,
)
from divbase_cli.user_auth import make_authenticated_request
from divbase_lib.api_schemas.queries import (
    BcftoolsQueryRequest,
    SampleMetadataQueryRequest,
    SampleMetadataQueryTaskResult,
)
from divbase_lib.api_schemas.task_history import TaskHistoryResult
from divbase_lib.exceptions import PolledTaskNotFinalError

logger = logging.getLogger(__name__)


METADATA_TSV_ARGUMENT = typer.Option(
    cli_settings.METADATA_TSV_NAME, help="Name of the sample metadata TSV file in the project's data store on DivBase."
)

BCFTOOLS_ARGUMENT = typer.Option(
    ...,
    help="""
        String consisting of the bcftools view command(s) to run. E.g. "view -r 21:15000000-25000000" or "view -s".
        The string cannot be empty; if you only want to subset on the selected samples, use: --command "view -s"
        """,
)

# Sample metadata and VCF queries both use the same core text, so it is defined up here.
TSV_FILTER_SYNTAX = (
    "String consisting of keys:values in the tsv file to filter on. "
    "The syntax is 'Key1:Value1,Value2;Key2:Value3,Value4', where the keys are the column header names in the tsv, "
    "and values are the column values. Multiple values for a key are separated by commas, and multiple keys are "
    "separated by semicolons. When multiple keys are provided, an intersect query will be performed. "
    "E.g. 'Area:West of Ireland,Northern Portugal;Sex:F'."
)

TSV_FILTER_HELP_TEXT_VCF = (
    "This option calculates the samples to filter the VCFs on based on a sample metadata query. "
    + TSV_FILTER_SYNTAX
    + "\n\nMutually exclusive with --samples, --samples-file, and --all-samples."
)

SAMPLE_SELECTION_HELP_PANEL = "Sample Selection (Required: Include Exactly One)"
VCF_QUERY_HELP_TEXT = (
    "Submit a VCF query to run on the DivBase API. "
    "A single, merged VCF file with the query results will be added to the project on success.\n\n"
    "Exactly one sample-selection mode is required: "
    "--tsv-filter | --samples | --samples-file | --all-samples."
)

GET_RESULTS_HELP_TEXT = (
    "Poll for completion of a query job and download (or print) the results."
    "Similar to running 'divbase-cli task-history id <TASK_ID>' but with the added benefit of polling for the"
    "terminal state of the job (SUCCESS/FAILED). Designed to be of particular use in scripts and other automated workflows."
)

query_app = typer.Typer(
    help="Run queries on the VCF files stored in the project's data store on DivBase. Queries are run on the DivBase API",
    no_args_is_help=True,
)


@query_app.command("tsv")
def sample_metadata_query(
    tsv_filter: str = typer.Argument(
        ...,
        help=TSV_FILTER_SYNTAX,
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
    """
    project_config = resolve_project(project_name=project)
    logged_in_url = ensure_logged_in(desired_url=project_config.divbase_url)

    request_data = SampleMetadataQueryRequest(tsv_filter=tsv_filter, metadata_tsv_name=metadata_tsv_name)

    response = make_authenticated_request(
        method="POST",
        divbase_base_url=logged_in_url,
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
            print(f"Sample ID: '{sample.sample_id}', Filename: '{sample.filename}'")

    print(f"The results for the query ([bright_blue]{results.query_message}[/bright_blue]):")

    unique_sample_ids = results.unique_sample_ids or []
    unique_filenames = results.unique_filenames or []
    print(f"Unique Sample IDs: {unique_sample_ids}")
    print(f"Unique filenames: {unique_filenames}\n")
    if not unique_sample_ids:
        print("[yellow]No samples match your query filters.[/yellow]\n")


@query_app.command("vcf", help=VCF_QUERY_HELP_TEXT)
def vcf_query(
    tsv_filter: str | None = typer.Option(
        None, help=TSV_FILTER_HELP_TEXT_VCF, rich_help_panel=SAMPLE_SELECTION_HELP_PANEL
    ),
    samples: str | None = typer.Option(
        None,
        help="Comma-separated list of sample IDs. Mutually exclusive with --tsv-filter, --samples-file, and --all-samples.",
        rich_help_panel=SAMPLE_SELECTION_HELP_PANEL,
    ),
    samples_file: Path | None = typer.Option(
        None,
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
        help=(
            "Path to a UTF-8 text file with one sample ID per line. Blank lines and lines starting with # are ignored. "
            "Mutually exclusive with --tsv-filter, --samples, and --all-samples."
        ),
        rich_help_panel=SAMPLE_SELECTION_HELP_PANEL,
    ),
    all_samples: bool = typer.Option(
        False,
        "--all-samples",
        help=(
            "Use all samples in the project for the query. "
            "Mutually exclusive with --tsv-filter, --samples, and --samples-file."
        ),
        rich_help_panel=SAMPLE_SELECTION_HELP_PANEL,
    ),
    command: str = BCFTOOLS_ARGUMENT,
    metadata_tsv_name: str = METADATA_TSV_ARGUMENT,
    project: str | None = PROJECT_NAME_OPTION,
) -> None:
    """
    Submit a VCF query to run on the DivBase API. A single, merged VCF file with the query results will be added to the project on success.
    """

    # TODO Error handling for subprocess calls.
    # TODO: handle the case empty results are returned from tsv_query()
    # TODO what if a job fails and the user wants to re-run it? do we store temp files?

    # Note! Pydantic model validator also enforces this on the API side just queries can be submitted directly to the endpoint.
    # This block here is to catch it on the CLI side with a more user-friendly error message before even making the API call.
    has_tsv_filter = tsv_filter is not None
    has_samples = samples is not None
    has_samples_file = samples_file is not None
    has_all_samples = all_samples
    selection_count = sum([has_tsv_filter, has_samples, has_samples_file, has_all_samples])
    if selection_count > 1:
        raise typer.BadParameter("Use only one of --tsv-filter, --samples, --samples-file, or --all-samples.")
    if selection_count == 0:
        raise typer.BadParameter(
            "Sample selection is required. Use one of --tsv-filter, --samples, --samples-file, or --all-samples."
        )

    normalized_samples, sample_input_warnings = _normalize_samples_input(samples=samples, samples_file=samples_file)
    if sample_input_warnings:
        print("[yellow]Warnings:[/yellow]")
        for warning in sample_input_warnings:
            print(f"  • {warning}")
        print()

    project_config = resolve_project(project_name=project)
    logged_in_url = ensure_logged_in(desired_url=project_config.divbase_url)

    request_data = BcftoolsQueryRequest(
        tsv_filter=tsv_filter,
        command=command,
        metadata_tsv_name=metadata_tsv_name,
        samples=normalized_samples,
        all_samples=all_samples,
    )

    response = make_authenticated_request(
        method="POST",
        divbase_base_url=logged_in_url,
        api_route=f"v1/query/vcf/projects/{project_config.name}",
        json=request_data.model_dump(),
    )

    task_id = response.json()
    print(
        f"Job submitted successfully with task id: {task_id}. To check the status of your job, use the command: divbase-cli task-history id {task_id}"
    )


@query_app.command("get-results", help=GET_RESULTS_HELP_TEXT)
def get_results_from_query_job_by_task_id(
    task_id: int = typer.Argument(..., help="Task ID of the query job to poll for results from."),
    project: str | None = PROJECT_NAME_OPTION,
) -> None:
    """
    Get results from a query job by its task ID by first polling for the final state of the task.
    """

    divbase_url = resolve_url_for_non_project_specific_commands()
    task_status = poll_task_until_final_state_reached(divbase_url=divbase_url, task_id=task_id)

    if task_status == "SUCCESS":
        print(f"Task {task_id} completed successfully.")
    else:
        print(f"Task {task_id} failed.")


@stamina.retry(
    on=retry_polling_until_final_or_retryable_api_errors,
    wait_initial=1.0,  # seconds. Some overhead time is required to init the Celery task even with idle workers, so wait a bit before first poll.
    wait_exp_base=2,  # exponential backoff factor
    wait_max=60.0,  # cap exponential backoff at once per 60 seconds if it reaches that far.
    timeout=60
    * 60
    * 12,  # To avoid infinite polling in case of unforeseen errors. This needs to include time in the queue and time to process the task
)
def poll_task_until_final_state_reached(divbase_url: str, task_id: int) -> TaskHistoryResult:
    """
    Poll for the final state (SUCCESS/FAILURE) of a celery task by task ID.

    Note that this is designed with two layers of retry_only_on_retryable_divbase_api_errors: the outer layer (this function), and the inner layer (make_authenticated_request; 3 attempts).
    This is to make this polling more resilient to temporary API errors for the intent of using the CLI command that calls this in scripts and other automated workflows.
    """

    FINAL_STATES = {"SUCCESS", "FAILURE"}

    response = make_authenticated_request(
        method="GET",
        divbase_base_url=divbase_url,
        api_route=f"v1/task-history/tasks/{task_id}",
    )
    # Endpoint returns 403 if the task is not found or the user does not have permission to view it. So item should not be empty.
    items = response.json()

    task_status = TaskHistoryResult(**items[0]).status  # for task ID lookups, only one entry is returned

    if task_status in FINAL_STATES:
        return task_status

    raise PolledTaskNotFinalError(f"Task {task_id} state is {task_status}")


def _normalize_samples_input(samples: str | None, samples_file: Path | None) -> tuple[list[str] | None, list[str]]:
    """
    Normalize sample selection inputs from CLI options.
    """
    if samples is not None:
        normalized = []
        for sample in samples.split(","):
            sample_clean = sample.strip()
            if sample_clean:
                normalized.append(sample_clean)

        if not normalized:
            raise typer.BadParameter("--samples must contain at least one non-empty sample ID.")
        return normalized, []

    if samples_file is not None:
        raw_lines = samples_file.read_text(encoding="utf-8").splitlines()
        normalized = []
        warnings = []
        delimiter_lines = []
        found_delimiters = set()
        disallowed_delimiters = (",", ";", "\t", "|")
        delimiter_display = {
            ",": "','",
            ";": "';'",
            "\t": "'tab'",
            "|": "'|'",
        }

        for line_number, raw_line in enumerate(raw_lines, start=1):
            stripped = raw_line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            normalized.append(stripped)
            line_delimiters = set()
            for delimiter in disallowed_delimiters:
                if delimiter in stripped:
                    line_delimiters.add(delimiter)

            if line_delimiters:
                delimiter_lines.append((line_number, stripped))
                found_delimiters.update(line_delimiters)

        if not normalized:
            raise typer.BadParameter(
                f"Samples file has no sample IDs after ignoring blank/comment lines: {samples_file}. Please ensure that it contains at least one sample ID."
            )

        if delimiter_lines:
            preview_entries = []  # Only preview up to the first 3 lines with delimiters to avoid overwhelming the user with a long list if there are many problematic lines.
            for line_number, line_value in delimiter_lines[:3]:
                preview_entries.append(f"line {line_number} ('{line_value}')")
            preview = ", ".join(preview_entries)

            previewed_line_count = min(len(delimiter_lines), 3)
            extra_line_count = len(delimiter_lines) - previewed_line_count
            extra_msg = f", +{extra_line_count} more line(s)" if extra_line_count > 0 else ""

            delimiter_names = []
            for delimiter in disallowed_delimiters:
                if delimiter in found_delimiters:
                    delimiter_names.append(delimiter_display[delimiter])
            delimiters_found_str = ", ".join(delimiter_names)

            raise typer.BadParameter(
                "Invalid --samples-file format: expected one sample ID per line with no delimiters. "
                f"Found delimiter(s) {delimiters_found_str} on {len(delimiter_lines)} line(s): \n{preview}{extra_msg}. "
                "\nPlease do not use delimiters (',' ';' '\\t' '|') in your samples file."
            )

        if len(normalized) == 1:
            warnings.append(
                "Only one sample ID was found in --samples-file after ignoring blank/comment lines. "
                "If this was not intended, verify that the file has one sample ID per line."
            )

        return normalized, warnings

    return None, []
