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
    # These are the states known by the worker. The state when a task is in the queue is handled by the broker and PENDING is typically used
    # To avoid having separate logic for the enqueued state, check task status against the worker state and return QUEUED to user's terminal.
    CELERY_STATES_EXCLUDING_PENDING = {"STARTED", "SUCCESS", "FAILURE", "RETRY", "REVOKED"}

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
        """Display the task history fetched from the results backend in a formatted table."""

        sorted_tasks = sorted(self.task_items.items(), key=lambda x: x[1].created_at or "", reverse=True)
        display_limit = self.display_limit or 10
        limited_tasks = sorted_tasks[:display_limit]

        table = self._create_task_history_table()

        for task_id, task in limited_tasks:
            raw_status = task.status
            if raw_status and raw_status.upper() in self.CELERY_STATES_EXCLUDING_PENDING:
                state = raw_status.upper()
            else:
                state = "QUEUING"

            colour = self.STATE_COLOURS.get(state, "white")
            state_with_colour = f"[{colour}]{state}[/{colour}]"

            submitter = task.submitter_email or "Unknown"
            result = self._format_result(task, state)

            table.add_row(
                submitter,
                task_id,
                state_with_colour,
                self._format_unix_timestamp(task.created_at),
                self._format_unix_timestamp(task.started_at),
                str(task.runtime if task.runtime is not None else "N/A"),
                result,
            )
        console = Console()
        console.print(table)

    def _format_unix_timestamp(self, timestamp):
        """
        The results backend task status returns timestamps as integers or floats.
        This function formats them into a human-readable string.
        """
        if timestamp is None:
            return "N/A"

        if isinstance(timestamp, datetime.datetime):
            dt = timestamp.replace(microsecond=0)
        elif isinstance(timestamp, (int, float)):
            dt = datetime.datetime.fromtimestamp(timestamp).replace(microsecond=0)
        else:
            return str(timestamp)

        if dt.tzinfo is not None:
            return dt.strftime("%Y-%m-%d %H:%M:%S %Z")
        else:
            local_dt = dt.astimezone()
            return local_dt.strftime("%Y-%m-%d %H:%M:%S %Z")

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
        table.add_column("Created at", style="yellow", width=19, overflow="fold")
        table.add_column("Started at", style="yellow", width=19, overflow="fold")
        table.add_column("Runtime (s)", style="blue", width=10, overflow="fold")
        table.add_column("Result", style="white", width=35, overflow="fold")
        return table

    def _format_result(self, task, state):
        """
        Format the result message based on the task state and type.
        """
        colour = self.STATE_COLOURS.get(state, "white")

        if state == "FAILURE":
            if isinstance(task.result, dict):
                error_msg = task.result.get("error")
                if not error_msg:
                    exc_message = task.result.get("exc_message")
                    if isinstance(exc_message, list) and exc_message:
                        error_msg = " ".join(str(msg) for msg in exc_message)
                    elif exc_message:
                        error_msg = str(exc_message)
                if not error_msg:
                    exc_type = task.result.get("exc_type")
                    error_msg = exc_type if exc_type else "Unknown error"
            else:
                error_msg = str(task.result) or "Unknown error"

            return f"[{colour}]{error_msg}[/{colour}]"

        if isinstance(task.result, BcftoolsQueryTaskResult):
            result_message = f"Output file ready for download: {task.result.output_file}"
            return f"[{colour}]{result_message}[/{colour}]"

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

        if isinstance(task.result, dict) and ("exc_type" in task.result or "error" in task.result):
            # Handle any remaining error dicts that weren't caught by FAILURE state check
            error_msg = task.result.get("error")
            if not error_msg:
                exc_message = task.result.get("exc_message")
                if isinstance(exc_message, list) and exc_message:
                    error_msg = " ".join(str(msg) for msg in exc_message)
                elif exc_message:
                    error_msg = str(exc_message)
                else:
                    error_msg = task.result.get("exc_type", "Unknown error")
            return f"[{colour}]{error_msg}[/{colour}]"

        result_message = str(task.result)
        return f"[{colour}]{result_message}[/{colour}]"
