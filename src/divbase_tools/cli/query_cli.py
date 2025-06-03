"""
Query subcommand for the divbase_tools CLI.
"""

import logging
from pathlib import Path

import typer
from rich import print

from divbase_tools.cli.user_config_cli import CONFIG_PATH_OPTION
from divbase_tools.cli.utils import resolve_bucket_name
from divbase_tools.cli.version_cli import BUCKET_NAME_OPTION
from divbase_tools.queries import dummy_pipe_query_command, pipe_query_command, tsv_query_command
from divbase_tools.services import download_files_command

logger = logging.getLogger(__name__)

query_app = typer.Typer(help="Query the metadata for the VCF files stored in the bucket.", no_args_is_help=True)


@query_app.command("tsv")
def tsv_query(
    file: Path = typer.Option(
        default=Path("./sample_metadata.tsv"),
        help="Path to the tsv metadata file.",
    ),
    filter: str = typer.Option(
        None,
        help="""
        String consisting of keys:values in the tsv file to filter on.
        The syntax is 'Key1:Value1,Value2;Key2:Value3,Value4', where the key
        are the column header names in the tsv, and values are the column values. 
        Multiple values for a key are separated by commas, and multiple keys are 
        separated by semicolons. When multple keys are provided, an intersect query 
        will be performed. E.g. 'Area:West of Ireland,Northern Portugal;Sex:F'.
        """,
    ),
    show_sample_results: bool = typer.Option(
        default=False,
        help="Print sample_ID and Filename results from the query.",
    ),
    bucket_name: str = BUCKET_NAME_OPTION,
    config_path: Path = CONFIG_PATH_OPTION,
) -> dict:
    # TODO it perhaps be useful to set the default download_dir in the config so that we can
    # look for files there? For now this code just uses file.parent as the download directory.

    # TODO handle when --filter = NONE

    # TODO handle when the name of the sample column is something other than Sample_ID

    if not file.exists():
        logger.info(f"No local copy of the tsv file found at: {file}. Checking bucket for file.")
        bucket_name = resolve_bucket_name(bucket_name, config_path)
        download_files_command(
            bucket_name=bucket_name,
            all_files=[file.name],
            download_dir=file.parent,
            bucket_version=None,
            config_path=config_path,
        )

    logger.info(f"Quering {file}\n")
    query_result, query_string = tsv_query_command(file=file, filter=filter)
    unique_sampleIDs = query_result["Sample_ID"].unique().tolist()
    unique_filenames = query_result["Filename"].unique().tolist()
    sample_and_filename_subset = query_result[["Sample_ID", "Filename"]]
    if show_sample_results:
        print("Name and file for each sample in query results:")
        print(f"{sample_and_filename_subset.to_string(index=False)}\n")

    print(f"The results for the query ([bright_blue]{query_string}[/bright_blue]):")
    print(f"Unique Sample IDs: {unique_sampleIDs}")
    print(f"Unique filenames: {unique_filenames}\n")

    return {
        "sample_and_filename_subset": sample_and_filename_subset,
        "sampleIDs": unique_sampleIDs,
        "filenames": unique_filenames,
    }


@query_app.command("bcftools-pipe")
def pipe_query(
    tsv_file: Path = typer.Option(
        default=Path("./sample_metadata.tsv"),
        help="Path to the tsv metadata file.",
    ),
    tsv_filter: str = typer.Option(
        None,
        help="""
        String consisting of keys:values in the tsv file to filter on.
        The syntax is 'Key1:Value1,Value2;Key2:Value3,Value4', where the key
        are the column header names in the tsv, and values are the column values. 
        Multiple values for a key are separated by commas, and multiple keys are 
        separated by semicolons. When multple keys are provided, an intersect query 
        will be performed. E.g. 'Area:West of Ireland,Northern Portugal;Sex:F'.
        """,
    ),
    command: str = typer.Option(
        None,
        help="""
        String consisting of the bcftools command to run on the files returned by the tsv query.
        """,
    ),
    bucket_name: str = BUCKET_NAME_OPTION,
    config_path: Path = CONFIG_PATH_OPTION,
    dry: bool = typer.Option(
        False,
        help="""Dry run that copies over a hard-coded query results file rather than generating it with bcftools""",
    ),
) -> None:
    # TODO Error handling for subprocess calls.
    # TODO: handle case empty results are returned from tsv_query()
    # TODO what if the user just want to run bcftools on existing files in the bucket, without a tsv file query first?
    # TODO what if a job fails and the user wants to re-run it? do we store temp files?
    if tsv_filter:
        unique_query_results = tsv_query(
            file=tsv_file,
            filter=tsv_filter,
            show_sample_results=False,
            bucket_name=bucket_name,
            config_path=config_path,
        )

    missing_files = []
    for filename in unique_query_results.get("filenames", []):
        file_path = Path.cwd() / filename
        if not file_path.exists():
            missing_files.append(filename)

    if missing_files:
        logger.warning(f"The following files were not found locally: {missing_files}")
        bucket_name = resolve_bucket_name(bucket_name, config_path)
        download_files_command(
            bucket_name=bucket_name,
            all_files=missing_files,
            download_dir=Path.cwd(),
            bucket_version=None,
            config_path=config_path,
        )

    if not dry:
        unique_query_results = (
            unique_query_results or {}
        )  # TODO handle this better, empty values will likely break downstream calls anyway
        pipe_query_command(command=command, bcftools_inputs=unique_query_results)
    else:
        dummy_pipe_query_command()
