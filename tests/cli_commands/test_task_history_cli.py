from contextlib import contextmanager
from time import sleep
from unittest.mock import patch

from typer.testing import CliRunner

from divbase_cli.display_task_history import TaskHistoryDisplayManager
from divbase_cli.divbase_cli import app
from tests.helpers.api_setup import TEST_USERS

runner = CliRunner()


@contextmanager
def capture_task_history_manager():
    """
    Context manager that captures the state of TaskHistoryDisplayManager instance so that
    assertions can be made on the actually content of the class and not the rich table output,
    which is sensitive to formatting errors (newlines, whitespaces, etc).
    """
    captured_manager = None

    def capture_display_manager(self):
        nonlocal captured_manager
        captured_manager = self
        return original_print_task_history(self)

    original_print_task_history = TaskHistoryDisplayManager.print_task_history

    with patch.object(TaskHistoryDisplayManager, "print_task_history", capture_display_manager):
        yield lambda: captured_manager

    assert captured_manager is not None, "No TaskHistoryDisplayManager was captured"


def test_edit_user_submits_and_sees_task_history(CONSTANTS, logged_in_edit_user_with_existing_config):
    """Integration test where edit user submits two tasks and sees them in their task history."""
    query_string = "Area:West of Ireland,Northern Portugal;Sex:F"
    project_name = CONSTANTS["QUERY_PROJECT"]
    result_submit = runner.invoke(app, f"dimensions update --project {project_name}")
    assert result_submit.exit_code == 0
    sleep(0.5)
    result_submit = runner.invoke(app, f"query tsv '{query_string}' --project {project_name}")
    assert result_submit.exit_code == 0

    with capture_task_history_manager() as get_manager:
        result_history = runner.invoke(app, "task-history user")
        assert result_history.exit_code == 0

        captured_manager = get_manager()

    assert captured_manager.command_context["user_name"] == TEST_USERS["edit user"]["email"]

    for task_id, task in captured_manager.task_items.items():
        print(task_id, task)
        assert task.kwargs.user_name == TEST_USERS["edit user"]["email"]
