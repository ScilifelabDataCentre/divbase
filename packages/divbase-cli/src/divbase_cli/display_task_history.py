import ast
import datetime
import logging

from rich.console import Console
from rich.table import Table

logger = logging.getLogger(__name__)


class TaskHistoryDisplayManager:
    """
    A class that manages displaying task history results to the user's terminal.
    """

    STATE_COLOURS = {
        "SUCCESS": "green",
        "FAILURE": "red",
        "PENDING": "yellow",
        "STARTED": "blue",
        "RETRY": "blue",
        "REVOKED": "magenta",
    }

    def __init__(self, task_items: dict, command_context: dict):
        self.task_items = task_items
        self.command_context = command_context

    def print_task_history(self) -> None:
        """Display the task history fetched from the Flower API in a formatted table."""

        sorted_tasks = sorted(self.task_items.items(), key=lambda x: x[1].state or "", reverse=True)
        display_limit = self.command_context.get("display_limit", 10)
        limited_tasks = sorted_tasks[:display_limit]

        table = self._create_task_history_table()

        for task_id, task in limited_tasks:
            state = task.state or "N/A"
            colour = self.STATE_COLOURS.get(state, "white")
            state_with_colour = f"[{colour}]{state}[/{colour}]"

            submitter = self._get_submitter_from_task_kwargs(task)
            result = self._format_result(task, state)

            table.add_row(
                submitter,
                task_id,
                state_with_colour,
                self._format_unix_timestamp(task.received),
                self._format_unix_timestamp(task.started),
                str(task.runtime if task.runtime is not None else "N/A"),
                result,
            )
        console = Console()
        console.print(table)

    def _format_unix_timestamp(self, timestamp):
        """
        The flower task status API returns timestamps as integers or floats.
        This function formats them into a human-readable string.
        """
        if isinstance(timestamp, (int, float)):
            dt = datetime.datetime.fromtimestamp(timestamp)
            local_timezone = datetime.datetime.now().astimezone().tzname()
            return f"{dt.strftime('%Y-%m-%d %H:%M:%S')} {local_timezone}"
        return str(timestamp)

    def _create_task_history_table(self):
        """
        Use the Rich library to initiate a table for displaying task history.
        """
        title_prefix = "DivBase Task History"
        if self.command_context.get("mode") == "user":
            title = f"{title_prefix} for User: {self.command_context.get('user_name', 'Unknown')}"
        elif self.command_context.get("mode") == "user_project":
            title = f"{title_prefix} for User: {self.command_context.get('user_name', 'Unknown')} and Project: {self.command_context.get('project_name', 'Unknown')}"
        elif self.command_context.get("mode") == "id":
            title = f"{title_prefix} for Task ID: {self.command_context.get('task_id', 'Unknown')}"
        elif self.command_context.get("mode") == "project":
            title = f"{title_prefix} for Project: {self.command_context.get('project_name', 'Unknown')}"
        else:
            title = title_prefix

        table = Table(title=title, show_lines=True)
        table.add_column("Submitting user", width=12, overflow="fold")
        table.add_column("Task ID", style="cyan")
        table.add_column("State", width=8)
        table.add_column("Received", style="yellow", width=19, overflow="fold")
        table.add_column("Started", style="yellow", width=19, overflow="fold")
        table.add_column("Runtime (s)", style="blue", width=10, overflow="fold")
        table.add_column("Result", style="white", width=35, overflow="fold")
        return table

    def _get_submitter_from_task_kwargs(self, task):
        """
        Extract submitter from task kwargs.
        """
        kwargs = task.kwargs or "{}"
        try:
            parsed_kwargs = ast.literal_eval(kwargs)
            return parsed_kwargs.get("user_name", "Unknown")
        except ValueError as e:
            logger.warning(f"Could not parse kwargs: {e}")
            return "Unknown"

    def _format_result(self, task, state):
        """
        Format the result message based on the task state.
        """
        colour = self.STATE_COLOURS.get(state, "white")
        if state == "FAILURE":
            exception_message = task.exception or "Unknown error"
            return f"[{colour}]{exception_message}[/{colour}]"
        else:
            result_message = str(task.result or "N/A")
            return f"[{colour}]{result_message}[/{colour}]"
