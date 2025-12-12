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

from divbase_cli.cli_commands.user_config_cli import CONFIG_FILE_OPTION
from divbase_cli.cli_commands.version_cli import PROJECT_NAME_OPTION
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
    config_file: Path = CONFIG_FILE_OPTION,
) -> None:
    """
    Query the tsv sidecar metadata file for the VCF files in the project's data store on DivBase.
    Returns the sample IDs and filenames that match the query.

    TODO: it perhaps be useful to set the default download_dir in the config so that we can
    look for files there? For now this code just uses file.parent as the download directory.
    TODO: handle when the name of the sample column is something other than Sample_ID
    """

    project_config = resolve_project(project_name=project, config_path=config_file)

    request_data = SampleMetadataQueryRequest(tsv_filter=filter, metadata_tsv_name=metadata_tsv_name)

    response = make_authenticated_request(
        method="POST",
        divbase_base_url=project_config.divbase_url,
        api_route=f"v1/query/sample-metadata/projects/{project_config.name}",
        json=request_data.model_dump(),
        timeout=20,  # This is longer than default (5), as api call response is query result, not a task-id.
    )

    results = SampleMetadataQueryTaskResult(**response.json())

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
    metadata_tsv_name: str = METADATA_TSV_ARGUMENT,
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
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
    project_config = resolve_project(project_name=project, config_path=config_file)

    request_data = BcftoolsQueryRequest(tsv_filter=tsv_filter, command=command, metadata_tsv_name=metadata_tsv_name)

    response = make_authenticated_request(
        method="POST",
        divbase_base_url=project_config.divbase_url,
        api_route=f"v1/query/bcftools-pipe/projects/{project_config.name}",
        json=request_data.model_dump(),
    )

    task_id = response.json()
    print(f"Job submitted successfully with task id: {task_id}")
