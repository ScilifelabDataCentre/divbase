"""
Task history subcommand for the DivBase CLI.

Submits a query for fetching the Celery task history for the user to the DivBase API.

"""

import logging

import typer

from divbase_cli.config_resolver import ensure_logged_in, resolve_project, resolve_url_for_non_project_specific_commands
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
    limit: int = typer.Option(10, help="Maximum number of tasks to display in the terminal. Sorted by recency."),
    project: str | None = typer.Option(
        None, help="Optional project name to filter the user's task history by project."
    ),
):
    """
    Check status of all tasks submitted by the user. Displays the latest 10 tasks by default, unless --limit is specified. Can be filtered by project name.
    """
    # TODO add option to sort ASC/DESC by task timestamp
    divbase_url = resolve_url_for_non_project_specific_commands()

    if project:
        task_history_response = make_authenticated_request(
            method="GET",
            divbase_base_url=divbase_url,
            api_route=f"v1/task-history/tasks/user/projects/{project}",
        )
    else:
        task_history_response = make_authenticated_request(
            method="GET",
            divbase_base_url=divbase_url,
            api_route="v1/task-history/tasks/user",
        )

    task_history_data = [TaskHistoryResult(**item) for item in task_history_response.json()]

    config = load_user_config()
    TaskHistoryDisplayManager(
        task_items=task_history_data,
        user_email=config.logged_in_email,
        project_name=project,
        mode="user_project" if project else "user",
        display_limit=limit,
    ).print_task_history()


@task_history_app.command("id")
def task_history_by_id(
    task_id: int | None = typer.Argument(..., help="Task ID to check the status of a specific query job."),
):
    """
    Check status of a specific task submitted by the user by its task ID.
    """
    divbase_url = resolve_url_for_non_project_specific_commands()

    task_history_response = make_authenticated_request(
        method="GET",
        divbase_base_url=divbase_url,
        api_route=f"v1/task-history/tasks/{task_id}",
    )
    task_history_data = [TaskHistoryResult(**item) for item in task_history_response.json()]

    TaskHistoryDisplayManager(
        task_items=task_history_data,
        user_email=None,
        project_name=None,
        mode="id",
    ).print_task_history()


@task_history_app.command("project")
def list_task_history_for_project(
    project: str = typer.Argument(
        None,
        help="Project name to check the task history for. Leave blank to use the default project set in your config.",
    ),
    limit: int = typer.Option(10, help="Maximum number of tasks to display in the terminal. Sorted by recency."),
):
    """
    Check status of all tasks submitted for a project. Requires a manager role in the project. Displays the latest 10 tasks by default, unless --limit is specified.
    """
    # TODO add option to sort ASC/DESC by task timestamp
    project_config = resolve_project(project_name=project)
    logged_in_url = ensure_logged_in(desired_url=project_config.divbase_url)

    task_history_response = make_authenticated_request(
        method="GET",
        divbase_base_url=logged_in_url,
        api_route=f"v1/task-history/projects/{project_config.name}",
    )

    task_history_data = [TaskHistoryResult(**item) for item in task_history_response.json()]

    TaskHistoryDisplayManager(
        task_items=task_history_data,
        user_email=None,
        project_name=project_config.name,
        mode="project",
        display_limit=limit,
    ).print_task_history()
