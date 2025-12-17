"""Unit tests for cron cleanup tasks."""

from unittest.mock import MagicMock, patch

import pytest

from divbase_api.worker.cron_tasks import cleanup_old_task_history_task, cleanup_stuck_tasks_task


@pytest.fixture
def mock_db_session():
    """Mock database session."""
    session = MagicMock()
    session.execute.return_value.rowcount = 5  # Simulates that 5 rows were deleted
    session.commit.return_value = None
    return session


def test_cleanup_old_task_history_deletes_old_entries(mock_db_session):
    """
    Test that cleanup_old_task_history_task deletes entries older than retention period.
    1. cleanup_old_task_history_task is run for a retention of 30 days
    2. The mock returns rowcount = 5, i.e. pretending 5 rows matched the delete statement

    """
    with patch("divbase_api.worker.cron_tasks.SyncSessionLocal") as mock_session_local:
        mock_session_local.return_value.__enter__.return_value = mock_db_session

        result = cleanup_old_task_history_task(retention_days=30)

        assert result["status"] == "completed"
        assert result["number_of_celery_meta_deleted"] == 5
        assert result["number_of_task_history_deleted"] == 5
        assert result["retention_days"] == 30
        assert mock_db_session.execute.call_count == 4
        assert mock_db_session.commit.called


def test_cleanup_stuck_tasks_deletes_stuck_pending_tasks(mock_db_session):
    """
    Test that cleanup_stuck_tasks removes tasks stuck in PENDING status.
    The mock db results simulates that there were 2 pending tasks and 0 started tasks
    that were found to be too old and were deleted.
    """
    mock_db_session.execute.return_value.fetchall.side_effect = [
        [("stuck-task-1",), ("stuck-task-2",)],  # PENDING tasks
        [],  # STARTED tasks
    ]

    with patch("divbase_api.worker.cron_tasks.SyncSessionLocal") as mock_session_local:
        mock_session_local.return_value.__enter__.return_value = mock_db_session

        result = cleanup_stuck_tasks_task(stuck_pending_hours=24, stuck_started_hours=48)

        assert result["status"] == "completed"
        assert result["number_of_stuck_pending_deleted"] == 2
        assert result["number_of_stuck_started_deleted"] == 0


def test_cleanup_stuck_tasks_deletes_stuck_started_tasks(mock_db_session):
    """Test that cleanup_stuck_tasks removes tasks stuck in STARTED status."""
    mock_db_session.execute.return_value.fetchall.side_effect = [
        [],  # PENDING tasks
        [("stuck-task-3",)],  # STARTED tasks
    ]

    with patch("divbase_api.worker.cron_tasks.SyncSessionLocal") as mock_session_local:
        mock_session_local.return_value.__enter__.return_value = mock_db_session

        result = cleanup_stuck_tasks_task(stuck_pending_hours=24, stuck_started_hours=48)

        assert result["status"] == "completed"
        assert result["number_of_stuck_pending_deleted"] == 0
        assert result["number_of_stuck_started_deleted"] == 1


def test_cleanup_tasks_handle_exceptions_gracefully():
    """Test that cleanup tasks handle database errors gracefully."""
    with patch("divbase_api.worker.cron_tasks.SyncSessionLocal") as mock_session_local:
        mock_session_local.return_value.__enter__.side_effect = Exception("Database connection error")

        with pytest.raises(Exception) as exc_info:
            cleanup_old_task_history_task(retention_days=30)

        assert "Database connection error" in str(exc_info.value)
