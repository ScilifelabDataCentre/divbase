"""
Query subcommand for the divbase_tools CLI.
"""

from pathlib import Path
import typer
from rich import print
import pandas as pd
from divbase_tools.queries import tsv_query_command, pipe_query_command

from divbase_tools.cli.user_config_cli import CONFIG_PATH_OPTION
from divbase_tools.cli.utils import resolve_bucket_name
from divbase_tools.cli.version_cli import BUCKET_NAME_OPTION
from divbase_tools.services import (
    download_files_command,
    list_files_command,
)
import logging
logger = logging.getLogger(__name__)

query_app = typer.Typer(
    help="Query the metadata for the VCF files stored in the bucket.", no_args_is_help=True
)

@query_app.command("tsv")
def tsv_query(
    file: Path = typer.Option(
        default=Path("./sample_metadata.tsv"),
        help="Path to the tsv metadata file.",
    ),
    filter: str = typer.Option(None, 
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
) -> None:
  
    # TODO it perhaps be useful to set the default download_dir in the config so that we can 
    # look for files there? For now this code just uses file.parent as the download directory.

    # TODO handle when --filter = NONE 

    # TODO handle when the name of the sample column is something other than Sample_ID
    
    if not file.exists():
        logger.info(f"No local copy of the tsv file found at: {file}. Checking bucket for file.")
        bucket_name = resolve_bucket_name(bucket_name, config_path)
        bucket_files = list_files_command(
            bucket_name=bucket_name,
            config_path=config_path,
        )
        if file.name in bucket_files:
            download_files_command(
                bucket_name=bucket_name,
                all_files=[file.name],
                download_dir=file.parent,
                bucket_version=None,
                config_path=config_path,
            )
        else:
            logger.error(f"{file} does not exist in the bucket. Exiting.")
            raise typer.Exit(1)

    logger.info(f"Quering {file}\n")
    query_result, query_string = tsv_query_command(file=file, filter=filter)
    unique_filenames = query_result["Filename"].unique().tolist()

    if show_sample_results:
        print(f"Name and file for each sample in query results:")  
        print(f"{query_result[["Sample_ID", "Filename"]].to_string(index=False)}\n")

    print(f"Unique filenames for query {query_string}: {unique_filenames}")

@query_app.command("bcftools-pipe")
def pipe_query(
    command: str = typer.Option(None, 
        help="""
        String consisting of the bcftools command to run on the files returned by the tsv query.
        """,
    ),
):

    #TODO Error handling for subprocess calles. Perhaps: subprocess.CalledProcessError
    pipe_query_command(command=command)