"""
E2E tests validating that task submission is blocked when the queue is closed,
and that admin users can still submit when the queue is closed.
"""

from datetime import datetime, timedelta, timezone

import pytest
from sqlalchemy import select
from typer.testing import CliRunner

from divbase_api.models.queue_status import QueueStatusDB
from divbase_cli.cli_exceptions import DivBaseAPIError
from divbase_cli.divbase_cli import app

runner = CliRunner()

PROJECT_NAME = "project1"
TSV_FILTER = "Area:West of Ireland"
BCFTOOLS_CMD = "view -r 21:15000000-25000000"
COMMANDS = [
    f"query vcf --tsv-filter '{TSV_FILTER}' --command '{BCFTOOLS_CMD}' --project {PROJECT_NAME}",
    f"query tsv '{TSV_FILTER}' --project {PROJECT_NAME}",
    f"dimensions update --project {PROJECT_NAME}",
]


@pytest.fixture
def queue_closed(db_session_sync):
    """
    Set the queue to closed before the test, restore it to open after.
    """
    result = db_session_sync.execute(select(QueueStatusDB).filter_by(id=1))
    queue_status = result.scalar_one()
    queue_status.is_closed = True
    queue_status.scheduled_start = None
    db_session_sync.commit()

    yield

    queue_status.is_closed = False
    queue_status.scheduled_start = None
    db_session_sync.commit()


@pytest.fixture
def queue_not_yet_closed(db_session_sync):
    """
    Set the queue to be closed in the future, jobs should still be able to be submitted.
    """
    result = db_session_sync.execute(select(QueueStatusDB).filter_by(id=1))
    queue_status = result.scalar_one()
    queue_status.is_closed = True
    queue_status.scheduled_start = datetime.now(tz=timezone.utc) + timedelta(days=1)
    db_session_sync.commit()

    yield

    queue_status.is_closed = False
    queue_status.scheduled_start = None
    db_session_sync.commit()


def test_regular_users_cannot_submit_tasks_when_queue_closed(
    CONSTANTS,
    logged_in_edit_user_with_existing_config,
    queue_closed,
):
    """Test that when the queue is closed, regular users cannot submit tasks"""
    for command in COMMANDS:
        result = runner.invoke(app, command)
        assert result.exit_code != 0
        assert isinstance(result.exception, DivBaseAPIError)
        assert "400" in str(result.exception)
        assert "queue_closed_error" in str(result.exception)


def test_admin_users_can_submit_tasks_when_queue_closed(
    CONSTANTS,
    logged_in_admin_with_existing_config,
    queue_closed,
):
    """Test that when the queue is closed, admin users can still submit tasks."""
    for command in COMMANDS:
        result = runner.invoke(app, command)
        # It's okay for the job to actually fail,
        # We just want the error to not be from the queue being closed
        # This check happens early enough, that a badly formulated query does not matter
        if result.exit_code != 0:
            assert "queue_closed_error" not in str(result.exception)


def test_regular_users_can_submit_tasks_when_queue_not_yet_closed(
    CONSTANTS,
    logged_in_edit_user_with_existing_config,
    queue_not_yet_closed,
):
    """Test that when the queue is not yet closed (just due to be closed), regular users can still submit tasks."""
    for command in COMMANDS:
        result = runner.invoke(app, command)
        # It's okay for the job to actually fail,
        # We just want the error to not be from the queue being closed
        # This check happens early enough, that a badly formulated query does not matter
        if result.exit_code != 0:
            assert "queue_closed_error" not in str(result.exception)
