import ast
import datetime
import logging
import os

import requests
from dotenv import load_dotenv
from rich.console import Console
from rich.table import Table

logger = logging.getLogger(__name__)

FLOWER_URL = "http://localhost:5555"


class TaskHistoryManager:
    """
    A class that manages interactions with the Flower API.
    It fetches, filters, and displays the task history based on the current user.
    """

    STATE_COLOURS = {
        "SUCCESS": "green",
        "FAILURE": "red",
        "PENDING": "yellow",
        "STARTED": "blue",
        "PROGRESS": "blue",
        "REVOKED": "magenta",
    }

    def __init__(self, url: str, flower_user: str = None, flower_password: str = None, divbase_user: str = None):
        self.url = url
        self.flower_API_user = flower_user
        self.flower_API_password = flower_password
        self.current_divbase_user = divbase_user

    def get_task_history(self, task_id: str = None, display_limit: int = 10) -> None:
        """
        Get the task history from the Flower API.
        """

        auth = (self.flower_API_user, self.flower_API_password)

        if task_id:
            request_url = f"{self.url}/api/task/info/{task_id}"
        else:
            api_limit = min(
                100, display_limit * 5
            )  # return more tasks than requested by the --limit arg to allow for downstream sorting
            request_url = f"{self.url}/api/tasks?limit={api_limit}"

        response = requests.get(request_url, auth=auth, timeout=3)

        if response.status_code == 200:
            self.print_task_history(response=response, task_id=task_id, display_limit=display_limit)
        else:
            print(f"Failed to fetch tasks: {response.status_code} - {response.text}")

    def print_task_history(self, response, task_id: str = None, display_limit: int = 10) -> None:
        """Display the task history fetched from the Flowwre APIin a formatted table."""
        tasks = response.json()

        task_items = []
        if task_id:
            task_items = [(task_id, tasks, self.parse_timestamp(tasks.get("started", 0)))]
        else:
            task_items = [(tid, data, self.parse_timestamp(data.get("started", 0))) for tid, data in tasks.items()]

        sorted_tasks = sorted(task_items, key=lambda x: x[2], reverse=True)
        limited_tasks = sorted_tasks[:display_limit]

        table = self.create_task_history_table()

        for id, task, _ in limited_tasks:
            state = task.get("state", "N/A")
            colour = self.STATE_COLOURS.get(state, "white")
            state_with_colour = f"[{colour}]{state}[/{colour}]"

            submitter = self.get_submitter_from_kwargs(task)

            result = self.format_result(task, state)

            if self.permission_to_view_task(submitter):
                table.add_row(
                    submitter,
                    id,
                    state_with_colour,
                    self.format_unix_timestamp(task.get("received", "N/A")),
                    self.format_unix_timestamp(task.get("started", "N/A")),
                    str(task.get("runtime", "N/A")),
                    result,
                )
        console = Console()
        console.print(table)

    def sort_tasks_by_runtime(self, tasks):
        """
        Sort tasks by their runtime in descending order.
        """
        return sorted(tasks, key=lambda x: x.get("runtime", 0), reverse=True)

    def parse_timestamp(self, timestamp):
        """
        Convert timestamp to float, if it is a numeric string.
        """
        if isinstance(timestamp, str) and timestamp.replace(".", "").isdigit():
            return float(timestamp)
        return timestamp

    def format_unix_timestamp(self, timestamp):
        """
        The flower task status API returns timestamps as integers or floats.
        This function formats them into a human-readable string.
        """
        if isinstance(timestamp, (int, float)):
            dt = datetime.datetime.fromtimestamp(timestamp)
            local_timezone = datetime.datetime.now().astimezone().tzname()
            return f"{dt.strftime('%Y-%m-%d %H:%M:%S')} {local_timezone}"
        return str(timestamp)

    def create_task_history_table(self):
        """
        Use the Rich library to initiate a table for displaying task history.
        """
        table = Table(title=f"DivBase Task Status for user: {self.current_divbase_user}", show_lines=True)
        table.add_column("Submitting user", width=12, overflow="fold")
        table.add_column("Task ID", style="cyan")
        table.add_column("State", width=8)
        table.add_column("Received", style="yellow", width=19, overflow="fold")
        table.add_column("Started", style="yellow", width=19, overflow="fold")
        table.add_column("Runtime (s)", style="blue", width=10, overflow="fold")
        table.add_column("Result", style="white", width=35, overflow="fold")
        return table

    def get_submitter_from_kwargs(self, task):
        """
        Extract submitter from task kwargs.
        """
        kwargs = task.get("kwargs", "{}")
        kwargs_dict = ast.literal_eval(kwargs)
        return kwargs_dict.get("submitter", "Unknown")

    def format_result(self, task, state):
        """
        Format the result message based on the task state.
        """
        if state == "FAILURE":
            exception_message = task.get("exception", "Unknown error")
            return f"[red]{exception_message}[/red]"
        else:
            result_message = str(task.get("result", "N/A"))
            return f"[green]{result_message}[/green]"

    def permission_to_view_task(self, submitter) -> bool:
        """
        Check if current user has permission to view this task
        """
        return self.current_divbase_user == submitter or self.current_divbase_user == "divbase_admin"


def dotenv_to_task_history_manager() -> TaskHistoryManager:
    """
    Create a TaskHistoryManager instance using the user's .env file.
    """
    load_dotenv()
    flower_user = os.environ.get("FLOWER_USER")
    flower_password = os.environ.get("FLOWER_PASSWORD")
    divbase_user = os.environ.get("DIVBASE_USER")

    if not flower_user:
        logger.warning("FLOWER_USER not provided or set in environment")
    if not flower_password:
        logger.warning("FLOWER_PASSWORD not provided or set in environment")
    if not divbase_user:
        logger.warning("DIVBASE_USER not set in environment")

    return TaskHistoryManager(
        url=FLOWER_URL,
        flower_user=flower_user,
        flower_password=flower_password,
        divbase_user=divbase_user,
    )
