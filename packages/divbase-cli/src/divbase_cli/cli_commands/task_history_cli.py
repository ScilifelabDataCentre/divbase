"""
Task history subcommand for the DivBase CLI.

Submits a query for fetching the Celery task history for the user to the DivBase API.

"""

import logging
from pathlib import Path

import typer

from divbase_cli.cli_commands.user_config_cli import CONFIG_FILE_OPTION
from divbase_cli.cli_exceptions import AuthenticationError
from divbase_cli.display_task_history import TaskHistoryDisplayManager
from divbase_cli.user_auth import make_authenticated_request
from divbase_cli.user_config import load_user_config
from divbase_lib.api_schemas.task_history import TaskHistoryResult

logger = logging.getLogger(__name__)


task_history_app = typer.Typer(
    help="Get the task history of query jobs submitted by the user to the DivBase API.",
    no_args_is_help=True,
)


@task_history_app.command("user")
def list_task_history_for_user(
    config_file: Path = CONFIG_FILE_OPTION,
    limit: int = typer.Option(10, help="Maximum number of tasks to display in the terminal. Sorted by recency."),
    project: str | None = typer.Option(
        None, help="Optional project name to filter the user's task history by project."
    ),
):
    """
    Check status of all tasks submitted by the user. Displays the latest 10 tasks by default, unless --limit is specified. Can be filtered by project name.
    """

    # TODO add option to sort ASC/DESC by task timestamp

    config = load_user_config(config_file)
    logged_in_url = config.logged_in_url

    if not logged_in_url:
        raise AuthenticationError("You are not logged in. Please log in with 'divbase-cli auth login [EMAIL]'.")

    if project:
        task_history_response = make_authenticated_request(
            method="GET",
            divbase_base_url=logged_in_url,
            api_route=f"v1/task-history/tasks/user/projects/{project}",
        )
    else:
        task_history_response = make_authenticated_request(
            method="GET",
            divbase_base_url=logged_in_url,
            api_route="v1/task-history/tasks/user",
        )

    task_history_data = [TaskHistoryResult(**item) for item in task_history_response.json()]

    TaskHistoryDisplayManager(
        task_items=task_history_data,
        user_name="TODO-GET-FROM-CONFIG-IN-THE-FUTURE",
        project_name=project,
        mode="user_project" if project else "user",
        display_limit=limit,
    ).print_task_history()


@task_history_app.command("id")
def task_history_by_id(
    task_id: int | None = typer.Argument(..., help="Task ID to check the status of a specific query job."),
    config_file: Path = CONFIG_FILE_OPTION,
):
    """
    Check status of a specific task submitted by the user by its task ID.
    """

    config = load_user_config(config_file)
    logged_in_url = config.logged_in_url

    if not logged_in_url:
        raise AuthenticationError("You are not logged in. Please log in with 'divbase-cli auth login [EMAIL]'.")

    task_history_response = make_authenticated_request(
        method="GET",
        divbase_base_url=logged_in_url,
        api_route=f"v1/task-history/tasks/{task_id}",
    )

    task_history_data = [TaskHistoryResult(**item) for item in task_history_response.json()]

    TaskHistoryDisplayManager(
        task_items=task_history_data,
        user_name=None,
        project_name=None,
        mode="id",
        display_limit=None,
    ).print_task_history()


@task_history_app.command("project")
def list_task_history_for_project(
    config_file: Path = CONFIG_FILE_OPTION,
    limit: int = typer.Option(10, help="Maximum number of tasks to display in the terminal. Sorted by recency."),
    project: str = typer.Argument(..., help="Project name to check the task history for."),
):
    """
    Check status of all tasks submitted for a project. Requires a manager role in the project. Displays the latest 10 tasks by default, unless --limit is specified.
    """

    # TODO add option to sort ASC/DESC by task timestamp
    # TODO use default project from config if not --project specified

    config = load_user_config(config_file)
    logged_in_url = config.logged_in_url

    if not logged_in_url:
        raise AuthenticationError("You are not logged in. Please log in with 'divbase-cli auth login [EMAIL]'.")

    task_history_response = make_authenticated_request(
        method="GET",
        divbase_base_url=logged_in_url,
        api_route=f"v1/task-history/projects/{project}",
    )

    task_history_data = [TaskHistoryResult(**item) for item in task_history_response.json()]

    TaskHistoryDisplayManager(
        task_items=task_history_data,
        user_name=None,
        project_name=project,
        mode="project",
        display_limit=limit,
    ).print_task_history()
