"""
Submit query style jobs to the DivBase API.
"""

from pathlib import Path

import httpx
import typer
from rich import print

from divbase_tools.cli_commands.query_cli import BCFTOOLS_ARGUEMENT, TSV_FILTER_ARGUEMENT
from divbase_tools.cli_commands.user_config_cli import CONFIG_FILE_OPTION
from divbase_tools.cli_commands.version_cli import BUCKET_NAME_OPTION
from divbase_tools.task_history import TaskHistoryManager
from divbase_tools.utils import resolve_bucket_name

DIVBASE_API_URL = "http://localhost:8000"


job_app = typer.Typer(help="Submit query jobs to the DivBase API via the DivBase CLI.", no_args_is_help=True)


@job_app.command("submit")
def submit_query_job(
    tsv_filter: str = TSV_FILTER_ARGUEMENT,
    command: str = BCFTOOLS_ARGUEMENT,
    bucket_name: str | None = BUCKET_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """Submit a query job (aka  bcftools-pipe) to the DivBase API."""
    bucket_name = resolve_bucket_name(bucket_name=bucket_name, config_path=config_file)

    params = {"tsv_filter": tsv_filter, "command": command, "bucket_name": bucket_name}
    response = httpx.post(f"{DIVBASE_API_URL}/jobs/", params=params)
    response.raise_for_status()

    task_id = response.json()
    print(f"Job submitted succsefully with task id: {task_id}")


@job_app.command("status")
def check_status(
    bucket_name: str | None = BUCKET_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """Check status of all jobs submitted by the user."""
    task_items = httpx.get(f"{DIVBASE_API_URL}/jobs/").json()
    task_history_manager = TaskHistoryManager(task_items=task_items, divbase_user="divbase_admin")
    task_history_manager.print_task_history()
