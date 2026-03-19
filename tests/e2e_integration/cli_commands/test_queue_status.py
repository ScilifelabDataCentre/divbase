"""
E2E tests validating that task submission is blocked when the queue is closed,
and that admin users can still submit when the queue is closed.
"""

import pytest
from sqlalchemy import select
from typer.testing import CliRunner

from divbase_api.models.queue_status import QueueStatus
from divbase_cli.cli_exceptions import DivBaseAPIError
from divbase_cli.divbase_cli import app

runner = CliRunner()


@pytest.fixture
def queue_closed(db_session_sync):
    """
    Set the queue to closed before the test, restore it to open after.
    """
    result = db_session_sync.execute(select(QueueStatus).filter_by(id=1))
    queue_status = result.scalar_one()
    queue_status.is_closed = True
    db_session_sync.commit()

    yield

    queue_status.is_closed = False
    db_session_sync.commit()


def test_regular_users_cannot_submit_tasks_when_queue_closed(
    CONSTANTS,
    logged_in_edit_user_with_existing_config,
    queue_closed,
):
    """Test that when the queue is closed, regular users cannot submit tasks"""

    project_name = CONSTANTS["QUERY_PROJECT"]
    tsv_filter = "Area:West of Ireland"
    arg_command = "view -s SAMPLES; view -r 21:15000000-25000000"

    commands = [
        f"query bcftools-pipe --tsv-filter '{tsv_filter}' --command '{arg_command}' --project {project_name}",
        f"query tsv '{tsv_filter}' --project {project_name}",
        f"dimensions update --project {project_name}",
    ]

    for command in commands:
        result = runner.invoke(app, command)
        assert result.exit_code != 0
        print(result.output)
        assert isinstance(result.exception, DivBaseAPIError)
        assert "400" in str(result.exception)
        assert "queue_closed_error" in str(result.exception)


def test_admin_users_can_submit_tasks_when_queue_closed(
    CONSTANTS,
    logged_in_admin_with_existing_config,
    queue_closed,
):
    """Test that when the queue is closed, admin users can still submit tasks."""
    project_name = CONSTANTS["QUERY_PROJECT"]
    tsv_filter = "Area:West of Ireland"
    arg_command = "view -s SAMPLES; view -r 21:15000000-25000000"

    commands = [
        f"query bcftools-pipe --tsv-filter '{tsv_filter}' --command '{arg_command}' --project {project_name}",
        f"query tsv Area:West of Ireland --project {project_name}",
        f"dimensions update --project {project_name}",
    ]

    for command in commands:
        result = runner.invoke(app, command)
        # It's okay for the job to actually fail,
        # We just want the error to not be from the queue being closed
        # This check happens early enough, that a badly formulated query does not matter
        if result.exit_code != 0:
            assert "queue_closed_error" not in str(result.exception)
