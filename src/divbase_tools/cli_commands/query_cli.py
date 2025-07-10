"""
Query subcommand for the divbase_tools CLI.

Submits queries (sample metadata and/or bcftools) to the DivBase API.

If sample metadata query:
    results are printed to the console.

If bcftools query:
    a task id is returned which can be used to check the status of the job.
    After task completed, a merged VCF file will be added to the bucket which can be downloaded by the user.


TODO:
- Ability to download results file given task id with the file cli?
-
"""

import logging
from pathlib import Path

import httpx
import typer
from rich import print

from divbase_tools.cli_commands.user_config_cli import CONFIG_FILE_OPTION
from divbase_tools.cli_commands.version_cli import BUCKET_NAME_OPTION
from divbase_tools.queries import SidecarQueryResult
from divbase_tools.task_history import TaskHistoryManager
from divbase_tools.utils import resolve_bucket_name

logger = logging.getLogger(__name__)


DIVBASE_API_URL = "http://localhost:8000"
DEFAULT_METADATA_TSV = Path("sample_metadata.tsv")

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
    help="Run queries on the VCF files stored in the bucket. Queries are run on the DivBase API", no_args_is_help=True
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
    bucket_name: str = BUCKET_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
) -> None:
    """
    Query the tsv sidecar metadata file for the VCF files stored in the bucket.
    Returns the sample IDs and filenames that match the query.

    TODO: it perhaps be useful to set the default download_dir in the config so that we can
    look for files there? For now this code just uses file.parent as the download directory.
    TODO: handle when the name of the sample column is something other than Sample_ID
    """
    bucket_name = resolve_bucket_name(bucket_name=bucket_name, config_path=config_file)

    params = {"tsv_filter": filter, "bucket_name": bucket_name}
    response = httpx.post(f"{DIVBASE_API_URL}/query/sample-metadata/", params=params)
    response.raise_for_status()

    results = SidecarQueryResult(**response.json())

    if show_sample_results:
        print("[bright_blue]Name and file for each sample in query results:[/bright_blue]")
        for sample in results.sample_and_filename_subset:
            print(f"Sample ID: '{sample['Sample_ID']}', Filename: '{sample['Filename']}'")

    print(f"The results for the query ([bright_blue]{results.query_message}[/bright_blue]):")
    print(f"Unique Sample IDs: {results.unique_sample_ids}")
    print(f"Unique filenames: {results.unique_filenames}\n")


@query_app.command("bcftools-pipe")
def pipe_query(
    tsv_filter: str = typer.Option(None, help=TSV_FILTER_HELP_TEXT),
 command: str = BCFTOOLS_ARGUMENT,
    bucket_name: str | None = BUCKET_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
) -> None:
    """
    Submit a query to run on the DivBase API. A single, merged VCF file will be added to the bucket on success.

    TODO Error handling for subprocess calls.
    TODO: handle case empty results are returned from tsv_query()
    TODO what if the user just want to run bcftools on existing files in the bucket, without a tsv file query first?
    TODO what if a job fails and the user wants to re-run it? do we store temp files?
    TODO be consistent about input argument and options. when are they optional, how is that indicated in docstring? etc.
    TODO consider handling the bcftools command whitelist checks also on the CLI level since the error messages are nicer looking?
    TODO consider moving downloading of missing files elsewhere, since this is now done before the celery task
    """
    bucket_name = resolve_bucket_name(bucket_name=bucket_name, config_path=config_file)

    params = {"tsv_filter": tsv_filter, "command": command, "bucket_name": bucket_name}
    response = httpx.post(f"{DIVBASE_API_URL}/query/bcftools-pipe/", params=params)
    response.raise_for_status()

    task_id = response.json()
    print(f"Job submitted succsefully with task id: {task_id}")


@query_app.command("task-status")
def check_status(
    task_id: str | None = typer.Argument(None, help="Optional task id to check the status of a specific query job."),
    config_file: Path = CONFIG_FILE_OPTION,
):
    """
    Check status of all query jobs submitted by the user.

    TODO - non bcftools query jobs should not show up here (but this should be handled by the API).
    TODO - Consider if representation of single task should be different from the list of tasks.
    """
    if task_id:
        task_items = httpx.get(f"{DIVBASE_API_URL}/query/{task_id}").json()
    else:
        task_items = httpx.get(f"{DIVBASE_API_URL}/query/").json()
    task_history_manager = TaskHistoryManager(task_items=task_items, divbase_user="divbase_admin")
    task_history_manager.print_task_history()
