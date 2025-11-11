import ast
import datetime
import logging

from rich.console import Console
from rich.table import Table

from divbase_lib.schemas.task_history import TaskHistoryResults

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

    # TODO use the command that user sent to the API as table header?

    def __init__(self, task_items: TaskHistoryResults):
        self.task_items = task_items

    def print_task_history(self, display_limit: int = 10) -> None:
        """Display the task history fetched from the Flower API in a formatted table."""

        sorted_tasks = sorted(self.task_items.tasks.items(), key=lambda x: x[1].state or "", reverse=True)
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
        table = Table(title="DivBase Task Status for TODO", show_lines=True)
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
