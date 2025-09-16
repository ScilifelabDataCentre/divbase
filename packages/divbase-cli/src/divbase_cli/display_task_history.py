import ast
import datetime
import logging

from rich.console import Console
from rich.table import Table

logger = logging.getLogger(__name__)


class TaskHistoryManager:
    """
    A class that manages interactions with the Flower API.
    It filters, and displays the task history based on the current user.
    """

    STATE_COLOURS = {
        "SUCCESS": "green",
        "FAILURE": "red",
        "PENDING": "yellow",
        "STARTED": "blue",
        "PROGRESS": "blue",
        "REVOKED": "magenta",
    }

    def __init__(self, task_items: list, divbase_user: str = None):
        self.task_items = task_items
        self.current_divbase_user = divbase_user

    def print_task_history(self, display_limit: int = 10) -> None:
        """Display the task history fetched from the Flowwre APIin a formatted table."""

        sorted_tasks = sorted(self.task_items, key=lambda x: x[2], reverse=True)
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
        # TODO, look into literal_eval.
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
