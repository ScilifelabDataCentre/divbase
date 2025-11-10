"""
Task history subcommand for the DivBase CLI.

Submits a query for fetching the Celery task history for the user to the DivBase API.

"""

import logging
from pathlib import Path

import typer

from divbase_cli.cli_commands.user_config_cli import CONFIG_FILE_OPTION
from divbase_cli.display_task_history import TaskHistoryManager
from divbase_cli.user_auth import make_authenticated_request
from divbase_cli.user_config import load_user_config
from divbase_lib.exceptions import AuthenticationError
from divbase_lib.queries import TaskHistoryResults

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

    params = {"limit": limit}

    if project:
        params["project_name"] = project
        task_history_response = make_authenticated_request(
            method="GET",
            divbase_base_url=logged_in_url,
            api_route="v1/task-history/list/user_and_project",
            params=params,
        )
    else:
        task_history_response = make_authenticated_request(
            method="GET",
            divbase_base_url=logged_in_url,
            api_route="v1/task-history/list/user",
            params=params,
        )

    task_history_data = TaskHistoryResults(**task_history_response.json())

    whoami_response = make_authenticated_request(
        method="GET",
        divbase_base_url=logged_in_url,
        api_route="v1/auth/whoami",
    )
    current_user_info = whoami_response.json()
    current_user_email = current_user_info["email"]
    is_admin = current_user_info.get("is_admin", False)

    # TODO this filters redundantly since the API should only return tasks for the current user unless the user is an admin
    task_history_manager = TaskHistoryManager(
        task_items=task_history_data, divbase_user=current_user_email, is_admin=is_admin
    )
    task_history_manager.print_task_history(display_limit=limit)


@task_history_app.command("id")
def task_history_by_id(
    task_id: str | None = typer.Argument(..., help="Task ID to check the status of a specific query job."),
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
        api_route=f"v1/task-history/list/{task_id}",
    )

    task_history_data = TaskHistoryResults(**task_history_response.json())

    whoami_response = make_authenticated_request(
        method="GET",
        divbase_base_url=logged_in_url,
        api_route="v1/auth/whoami",
    )
    current_user_info = whoami_response.json()
    current_user_email = current_user_info["email"]
    is_admin = current_user_info.get("is_admin", False)

    task_history_manager = TaskHistoryManager(
        task_items=task_history_data, divbase_user=current_user_email, is_admin=is_admin
    )
    task_history_manager.print_task_history()


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

    params = {"limit": limit}
    if project:
        params["project_name"] = project

    task_history_response = make_authenticated_request(
        method="GET",
        divbase_base_url=logged_in_url,
        api_route="v1/task-history/project",
        params=params,
    )

    task_history_data = TaskHistoryResults(**task_history_response.json())

    whoami_response = make_authenticated_request(
        method="GET",
        divbase_base_url=logged_in_url,
        api_route="v1/auth/whoami",
    )
    current_user_info = whoami_response.json()
    current_user_email = current_user_info["email"]
    is_admin = current_user_info.get("is_admin", False)

    task_history_manager = TaskHistoryManager(
        task_items=task_history_data, divbase_user=current_user_email, is_admin=is_admin
    )
    task_history_manager.print_task_history(display_limit=limit)
