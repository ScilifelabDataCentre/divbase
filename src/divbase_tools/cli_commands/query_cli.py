"""
Query subcommand for the divbase_tools CLI.
"""

import logging
import os
from pathlib import Path

import typer
from dotenv import load_dotenv
from rich import print

from divbase_tools.cli_commands.user_config_cli import CONFIG_FILE_OPTION
from divbase_tools.cli_commands.version_cli import BUCKET_NAME_OPTION
from divbase_tools.queries import SidecarQueryManager, fetch_query_files_from_bucket
from divbase_tools.task_history import TaskHistoryManager, get_task_history
from divbase_tools.tasks import bcftools_pipe_task
from divbase_tools.utils import resolve_bucket_name

logger = logging.getLogger(__name__)

query_app = typer.Typer(help="Query the VCF files stored in the bucket.", no_args_is_help=True)


TSV_FILTER_ARGUEMENT = typer.Option(
    None,
    help="""
    String consisting of keys:values in the tsv file to filter on.
    The syntax is 'Key1:Value1,Value2;Key2:Value3,Value4', where the key
    are the column header names in the tsv, and values are the column values. 
    Multiple values for a key are separated by commas, and multiple keys are 
    separated by semicolons. When multple keys are provided, an intersect query 
    will be performed. E.g. 'Area:West of Ireland,Northern Portugal;Sex:F'.
    """,
)

BCFTOOLS_ARGUEMENT = typer.Option(
    ...,
    help="""
        String consisting of the bcftools command to run on the files returned by the tsv query.
        """,
)


@query_app.command("tsv")
def tsv_query(
    file: Path = typer.Option(
        default=Path("./sample_metadata.tsv"),
        help="Path to the tsv metadata file.",
    ),
    filter: str = typer.Argument(
        ...,
        help="""
        String consisting of keys:values in the tsv file to filter on.
        The syntax is 'Key1:Value1,Value2;Key2:Value3,Value4', where the key
        are the column header names in the tsv, and values are the column values. 
        Multiple values for a key are separated by commas, and multiple keys are 
        separated by semicolons. When multple keys are provided, an intersect query 
        will be performed. E.g. 'Area:West of Ireland,Northern Portugal;Sex:F'.
        Use '' (empty string in quotes) to return ALL records.
        """,
    ),
    show_sample_results: bool = typer.Option(
        default=False,
        help="Print sample_ID and Filename results from the query.",
    ),
    bucket_name: str = BUCKET_NAME_OPTION,
    config_path: Path = CONFIG_FILE_OPTION,
) -> dict:
    """Query the tsv sidecar metadata file for the VCF files stored in the bucket. Returns the sample IDs and filenames that match the query."""
    # TODO it perhaps be useful to set the default download_dir in the config so that we can
    # look for files there? For now this code just uses file.parent as the download directory.

    # TODO handle when the name of the sample column is something other than Sample_ID

    if not file.exists():
        logger.info(f"No local copy of the tsv file found at: {file}. Checking bucket for file.")
        fetch_query_files_from_bucket(
            bucket_name=bucket_name, config_path=config_path, files=[file.name], download_dir=None, bucket_version=None
        )

    logger.info(f"Querying {file}\n")
    sidecar_manager = SidecarQueryManager(file=file).run_query(filter_string=filter)
    query_result = sidecar_manager.query_result
    query_message = sidecar_manager.query_message

    unique_sampleIDs = sidecar_manager.get_unique_values("Sample_ID")
    unique_filenames = sidecar_manager.get_unique_values("Filename")
    sample_and_filename_subset = query_result[["Sample_ID", "Filename"]]
    serialized_samples = sample_and_filename_subset.to_dict(orient="records")

    if show_sample_results:
        print("Name and file for each sample in query results:")
        print(f"{sample_and_filename_subset.to_string(index=False)}\n")

    print(f"The results for the query ([bright_blue]{query_message}[/bright_blue]):")
    print(f"Unique Sample IDs: {unique_sampleIDs}")
    print(f"Unique filenames: {unique_filenames}\n")

    return {
        "sample_and_filename_subset": serialized_samples,
        "sampleIDs": unique_sampleIDs,
        "filenames": unique_filenames,
    }


@query_app.command("bcftools-pipe")
def pipe_query(
    tsv_file: Path = typer.Option(
        default=Path("./sample_metadata.tsv"),
        help="Path to the tsv metadata file.",
    ),
    tsv_filter: str = TSV_FILTER_ARGUEMENT,
    command: str = BCFTOOLS_ARGUEMENT,
    bucket_name: str = BUCKET_NAME_OPTION,
    config_path: Path = CONFIG_FILE_OPTION,
    run_async: bool = typer.Option(False, "--async", help="Run as async job using Celery"),
) -> None:
    """
    Run bcftools commands on the files returned by an optional tsv query. Returns a merged VCF file.
    """
    # TODO Error handling for subprocess calls.
    # TODO: handle case empty results are returned from tsv_query()
    # TODO what if the user just want to run bcftools on existing files in the bucket, without a tsv file query first?
    # TODO what if a job fails and the user wants to re-run it? do we store temp files?
    # TODO be consistent about input argument and options. when are they optional, how is that indicated in docstring? etc.
    # TODO consider handling the bcftools command whitelist checks also on the CLI level since the error messages are nicer looking?
    # TODO consider moving downloading of missing files elsewhere, since this is now done before the celery task
    bucket_name = resolve_bucket_name(bucket_name=bucket_name, config_path=config_path)

    if not command or command.strip() == "" or command.strip() == ";":
        logger.error("Empty command provided. Please specify at least one valid bcftools command.")
        raise typer.Exit(code=1)

    filter = tsv_filter if tsv_filter else ""

    unique_query_results = tsv_query(
        file=tsv_file,
        filter=filter,
        show_sample_results=False,
        bucket_name=bucket_name,
        config_path=config_path,
    )
    if not unique_query_results.get("filenames"):
        logger.error("No files found matching your query criteria. Cannot proceed with bcftools commands.")
        raise typer.Exit(code=1)

    if not unique_query_results.get("sampleIDs"):
        logger.error("No samples found matching your query criteria. Cannot proceed with bcftools commands.")
        raise typer.Exit(code=1)

    if not tsv_filter:
        logger.info(f"No filter provided - using all {len(unique_query_results['sampleIDs'])} samples from {tsv_file}")

    missing_files = []
    for filename in unique_query_results.get("filenames", []):
        file_path = Path.cwd() / filename
        if not file_path.exists():
            missing_files.append(filename)
    if missing_files:
        logger.warning(f"The following files were not found locally: {missing_files}")
        fetch_query_files_from_bucket(
            bucket_name=bucket_name,
            config_path=config_path,
            files=missing_files,
            download_dir=None,
            bucket_version=None,
        )

    bcftools_inputs = unique_query_results

    load_dotenv()
    current_divbase_user = os.environ.get("DIVBASE_USER")
    if not current_divbase_user:
        logger.error("DIVBASE_USER environment variable not set")

    task_kwargs = {
        "command": command,
        "bcftools_inputs": bcftools_inputs,
        "submitter": current_divbase_user,
        "bucket_name": bucket_name,
    }

    if run_async:
        result = bcftools_pipe_task.apply_async(kwargs=task_kwargs)
        print(f"Job submitted with task ID: {result.id}")
        return result.id
    else:
        result = bcftools_pipe_task.apply(kwargs=task_kwargs)
        return result.id


@query_app.command("task-status")
def celery_task_status(
    task_id: str = typer.Option(
        None,
        help="Optional: task ID for the Celery task to get the status of.",
    ),
    limit: int = typer.Option(10, help="Optional: number of tasks to show (when not filtering by ID)"),
) -> None:
    """
    Check the query task history for the current user, either by task ID or by showing the last N tasks.
    """
    task_items = get_task_history(task_id=task_id, display_limit=limit)
    task_history_manager = TaskHistoryManager(task_items=task_items, divbase_user=os.environ.get("DIVBASE_USER"))
    task_history_manager.print_task_history()
