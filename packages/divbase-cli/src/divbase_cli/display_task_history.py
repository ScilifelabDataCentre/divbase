import datetime
import logging

from rich.console import Console
from rich.table import Table

from divbase_lib.api_schemas.task_history import (
    BcftoolsQueryTaskResult,
    DimensionUpdateTaskResult,
    SampleMetadataQueryTaskResult,
)

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

    def __init__(
        self,
        task_items: dict,
        user_name: str | None,
        project_name: str | None,
        task_id: str | None,
        mode: str,
        display_limit: int,
    ):
        self.task_items = task_items
        self.user_name = user_name
        self.project_name = project_name
        self.task_id = task_id
        self.mode = mode
        self.display_limit = display_limit

    def print_task_history(self) -> None:
        """Display the task history fetched from the Flower API in a formatted table."""

        sorted_tasks = sorted(self.task_items.items(), key=lambda x: x[1].state or "", reverse=True)
        display_limit = self.display_limit or 10
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
        if self.mode == "user":
            title = f"{title_prefix} for User: {self.user_name or 'Unknown'}"
        elif self.mode == "user_project":
            title = (
                f"{title_prefix} for User: {self.user_name or 'Unknown'} and Project: {self.project_name or 'Unknown'}"
            )
        elif self.mode == "id":
            title = f"{title_prefix} for Task ID: {self.task_id or 'Unknown'}"
        elif self.mode == "project":
            title = f"{title_prefix} for Project: {self.project_name or 'Unknown'}"
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
        # TODO decide if there are better ways of getting the submitting user than from the task kwargs. a lookup in the task_history db table would maybe make more sense?
        kwargs = task.kwargs
        if kwargs is None:
            return "Unknown"
        user_name = getattr(kwargs, "user_name", None)
        if user_name:
            return user_name
        if isinstance(kwargs, dict):
            return kwargs.get("user_name", "Unknown")
        return "Unknown"

    def _format_result(self, task, state):
        """
        Format the result message based on the task state and type.
        """
        colour = self.STATE_COLOURS.get(state, "white")
        if state == "FAILURE":
            exception_message = task.exception or "Unknown error"
            return f"[{colour}]{exception_message}[/{colour}]"

        if isinstance(task.result, BcftoolsQueryTaskResult):
            result_message = f"Output file ready for download: {task.result.output_file}"
            return f"[{colour}]{result_message}[/{colour}]"
            # TODO this should maybe say which project?

        elif isinstance(task.result, SampleMetadataQueryTaskResult):
            result_message = (
                f"Unique sample IDs:\n  {task.result.unique_sample_ids}\n"
                f"VCF files containing the sample IDs:\n  {task.result.unique_filenames}\n"
                f"Sample metadata query:\n  {task.result.query_message}"
            )
            return f"[{colour}]{result_message}[/{colour}]"

        elif isinstance(task.result, DimensionUpdateTaskResult):
            result_message = (
                f"VCF file dimensions index added or updated:\n  {task.result.VCF_files_added}\n"
                f"VCF files skipped by this job (previous DivBase-generated result VCFs):\n  {task.result.VCF_files_skipped}\n"
                f"VCF files that have been deleted from the project and now are dropped from the index:\n  {task.result.VCF_files_deleted}"
            )
            return f"[{colour}]{result_message}[/{colour}]"

        result_message = str(task.result)
        return f"[{colour}]{result_message}[/{colour}]"
