"""
Integration tests for cron cleanup tasks.

Tests the actual cleanup logic by creating database entries with backdated timestamps
and verifying that the cleanup tasks correctly delete old entries.
"""

import time
from datetime import datetime, timedelta, timezone

import pytest
from sqlalchemy import select, text

from divbase_api.models.task_history import CeleryTaskMeta, TaskHistoryDB, TaskStatus
from divbase_api.worker.cron_tasks import cleanup_old_task_history_task, cleanup_stuck_tasks_task


@pytest.fixture
def create_task_entry(db_session_sync, project_map):
    """
    Fixture to create task history entries with backdated timestamps.
    Returns a function that creates entries at specified ages.
    """

    def _create_entry(
        time_old: int,
        time_unit: str = "days",
        status: TaskStatus = TaskStatus.SUCCESS,
        project_id: int | None = "default",
        user_id: int = 1,
    ) -> str:
        """
        Create a task entry in both TaskHistoryDB and CeleryTaskMeta with backdated timestamps.

        Since CeleryTaskMeta is managed by Celery, there were issues with autoincrementing ID when
        trying to generate db entries with SQLalchemy. Thus, raw SQL is used to insert into
        CeleryTaskMeta with explicit IDs.
        """
        if project_id == "default":
            project_id = list(project_map.values())[0]

        # There are cron tasks that use days, and hours, so support both units
        if time_unit == "days":
            delta_time = datetime.now(timezone.utc) - timedelta(days=time_old)
        elif time_unit == "hours":
            delta_time = datetime.now(timezone.utc) - timedelta(hours=time_old)
        else:
            raise ValueError(f"Invalid time_unit: {time_unit}. Must be 'days' or 'hours'")

        timestamp = int(time.time() * 1000000)

        prefix = "system-task" if project_id is None else "task"
        task_id = f"{prefix}-{time_old}{time_unit[0]}-{status.value}-{timestamp}"

        # Determine timestamps based on status
        # For PENDING: only created_at is set
        # For STARTED: created_at and started_at are set
        # For SUCCESS/FAILURE: all three timestamps are set
        created_at = delta_time
        started_at = delta_time if status in [TaskStatus.STARTED, TaskStatus.SUCCESS, TaskStatus.FAILURE] else None
        completed_at = delta_time if status in [TaskStatus.SUCCESS, TaskStatus.FAILURE] else None

        task_history = TaskHistoryDB(
            task_id=task_id,
            user_id=user_id,
            project_id=project_id,
            status=status,
            created_at=created_at,
            started_at=started_at,
            completed_at=completed_at,
        )
        db_session_sync.add(task_history)
        db_session_sync.commit()

        # Create CeleryTaskMeta entry using raw SQL with explicit id
        max_id_result = db_session_sync.execute(text("SELECT COALESCE(MAX(id), 0) + 1 FROM celery_taskmeta")).scalar()
        db_session_sync.execute(
            text("""
                INSERT INTO celery_taskmeta 
                (id, task_id, status, result, date_done, traceback, name, args, kwargs, worker, retries, queue)
                VALUES 
                (:id, :task_id, :status, :result, :date_done, :traceback, :name, :args, :kwargs, :worker, :retries, :queue)
            """),
            {
                "id": max_id_result,
                "task_id": task_id,
                "status": status.value.upper(),
                "result": b'{"status": "completed"}' if status in [TaskStatus.SUCCESS, TaskStatus.FAILURE] else None,
                "date_done": completed_at,
                "traceback": None,
                "name": None,
                "args": None,
                "kwargs": None,
                "worker": None,
                "retries": None,
                "queue": None,
            },
        )
        db_session_sync.commit()

        return task_id

    return _create_entry


def test_cleanup_old_task_history_deletes_entries_older_than_threshold(
    db_session_sync,
    create_task_entry,
):
    """
    Test that cleanup_old_task_history_task deletes entries older than retention period
    and keeps recent entries.
    """
    retention_days = 30

    # Create old entries (should be deleted)
    old_task_id_35days = create_task_entry(time_old=35, status=TaskStatus.SUCCESS)
    old_task_id_50days = create_task_entry(time_old=50, status=TaskStatus.FAILURE)

    # Create recent entries (should be kept)
    recent_task_id_10days = create_task_entry(time_old=10, status=TaskStatus.SUCCESS)
    recent_task_id_25days = create_task_entry(time_old=25, status=TaskStatus.SUCCESS)

    result = cleanup_old_task_history_task(retention_days=retention_days)

    assert result["status"] == "completed"
    assert result["number_of_celery_meta_deleted"] == 2
    assert result["number_of_task_history_deleted"] == 2
    assert result["retention_days"] == retention_days

    for task_id in [old_task_id_35days, old_task_id_50days]:
        task_history = db_session_sync.execute(
            select(TaskHistoryDB).where(TaskHistoryDB.task_id == task_id)
        ).scalar_one_or_none()
        assert task_history is None, f"Task {task_id} should have been deleted from TaskHistoryDB"

        celery_meta = db_session_sync.execute(
            select(CeleryTaskMeta).where(CeleryTaskMeta.task_id == task_id)
        ).scalar_one_or_none()
        assert celery_meta is None, f"Task {task_id} should have been deleted from CeleryTaskMeta"

    for task_id in [recent_task_id_10days, recent_task_id_25days]:
        task_history = db_session_sync.execute(
            select(TaskHistoryDB).where(TaskHistoryDB.task_id == task_id)
        ).scalar_one_or_none()
        assert task_history is not None, f"Task {task_id} should have been kept in TaskHistoryDB"

        celery_meta = db_session_sync.execute(
            select(CeleryTaskMeta).where(CeleryTaskMeta.task_id == task_id)
        ).scalar_one_or_none()
        assert celery_meta is not None, f"Task {task_id} should have been kept in CeleryTaskMeta"


def test_cleanup_old_task_history_with_no_old_entries(db_session_sync, create_task_entry):
    """
    Test that cleanup task handles the case where there are no old entries gracefully.
    """
    retention_days = 30

    create_task_entry(time_old=5, status=TaskStatus.SUCCESS)
    create_task_entry(time_old=15, status=TaskStatus.SUCCESS)

    result = cleanup_old_task_history_task(retention_days=retention_days)

    assert result["status"] == "completed"
    assert result["number_of_celery_meta_deleted"] == 0
    assert result["number_of_task_history_deleted"] == 0


def test_cleanup_old_task_history_with_system_tasks(db_session_sync, create_task_entry):
    """
    Test that cleanup works correctly for system tasks (project_id=NULL).
    """
    retention_days = 30

    system_task_id = create_task_entry(
        time_old=35,
        status=TaskStatus.SUCCESS,
        project_id=None,
        user_id=2,
    )

    result = cleanup_old_task_history_task(retention_days=retention_days)

    assert result["status"] == "completed"
    assert result["number_of_task_history_deleted"] >= 1
    assert result["number_of_celery_meta_deleted"] >= 1

    task_history = db_session_sync.execute(
        select(TaskHistoryDB).where(TaskHistoryDB.task_id == system_task_id)
    ).scalar_one_or_none()
    assert task_history is None

    celery_meta = db_session_sync.execute(
        select(CeleryTaskMeta).where(CeleryTaskMeta.task_id == system_task_id)
    ).scalar_one_or_none()
    assert celery_meta is None


def test_cleanup_stuck_tasks_removes_pending_and_started(db_session_sync, create_task_entry):
    """
    Test that cleanup_stuck_tasks removes PENDING and STARTED tasks stuck longer than their respective thresholds.
    """
    stuck_pending_hours = 24
    stuck_started_hours = 48

    # Create stuck tasks that should be deleted
    stuck_pending_30h = create_task_entry(time_old=30, time_unit="hours", status=TaskStatus.PENDING)
    stuck_pending_50h = create_task_entry(time_old=50, time_unit="hours", status=TaskStatus.PENDING)
    stuck_started_50h = create_task_entry(time_old=50, time_unit="hours", status=TaskStatus.STARTED)
    stuck_started_72h = create_task_entry(time_old=72, time_unit="hours", status=TaskStatus.STARTED)

    # Create recent tasks that should be kept
    recent_pending_10h = create_task_entry(time_old=10, time_unit="hours", status=TaskStatus.PENDING)
    recent_started_30h = create_task_entry(time_old=30, time_unit="hours", status=TaskStatus.STARTED)

    result = cleanup_stuck_tasks_task(
        stuck_pending_hours=stuck_pending_hours,
        stuck_started_hours=stuck_started_hours,
    )

    assert result["status"] == "completed"
    assert result["number_of_stuck_pending_deleted"] == 2
    assert result["number_of_stuck_started_deleted"] == 2

    for task_id in [stuck_pending_30h, stuck_pending_50h, stuck_started_50h, stuck_started_72h]:
        task_history = db_session_sync.execute(
            select(TaskHistoryDB).where(TaskHistoryDB.task_id == task_id)
        ).scalar_one_or_none()
        assert task_history is None, f"Stuck task {task_id} should have been deleted from TaskHistoryDB"

        celery_meta = db_session_sync.execute(
            select(CeleryTaskMeta).where(CeleryTaskMeta.task_id == task_id)
        ).scalar_one_or_none()
        assert celery_meta is None, f"Stuck task {task_id} should have been deleted from CeleryTaskMeta"

    for task_id in [recent_pending_10h, recent_started_30h]:
        task_history = db_session_sync.execute(
            select(TaskHistoryDB).where(TaskHistoryDB.task_id == task_id)
        ).scalar_one_or_none()
        assert task_history is not None, f"Task {task_id} should have been kept in TaskHistoryDB"

        celery_meta = db_session_sync.execute(
            select(CeleryTaskMeta).where(CeleryTaskMeta.task_id == task_id)
        ).scalar_one_or_none()
        assert celery_meta is not None, f"Task {task_id} should have been kept in CeleryTaskMeta"


def test_cleanup_stuck_tasks_preserves_completed_tasks(db_session_sync, create_task_entry):
    """
    Test that stuck task cleanup only targets PENDING/STARTED, not SUCCESS/FAILURE.
    """
    stuck_pending_hours = 1
    stuck_started_hours = 1

    # Create old completed tasks (should NOT be deleted by stuck task cleanup. They are handled by old task history cleanup which is not tested here)
    old_success_id = create_task_entry(time_old=100, time_unit="hours", status=TaskStatus.SUCCESS)
    old_failure_id = create_task_entry(time_old=100, time_unit="hours", status=TaskStatus.FAILURE)

    result = cleanup_stuck_tasks_task(
        stuck_pending_hours=stuck_pending_hours,
        stuck_started_hours=stuck_started_hours,
    )

    assert result["status"] == "completed"

    for task_id in [old_success_id, old_failure_id]:
        task_history = db_session_sync.execute(
            select(TaskHistoryDB).where(TaskHistoryDB.task_id == task_id)
        ).scalar_one_or_none()
        assert task_history is not None, f"Completed task {task_id} should not have been deleted from TaskHistoryDB"

        celery_meta = db_session_sync.execute(
            select(CeleryTaskMeta).where(CeleryTaskMeta.task_id == task_id)
        ).scalar_one_or_none()
        assert celery_meta is not None, f"Completed task {task_id} should not have been deleted from CeleryTaskMeta"


def test_cleanup_stuck_tasks_deletes_from_both_tables(db_session_sync, create_task_entry):
    """
    Test that cleanup_stuck_tasks deletes from both TaskHistoryDB and CeleryTaskMeta.
    """
    stuck_pending_hours = 24
    stuck_started_hours = 48

    stuck_task_id = create_task_entry(time_old=30, time_unit="hours", status=TaskStatus.PENDING)

    task_history_before = db_session_sync.execute(
        select(TaskHistoryDB).where(TaskHistoryDB.task_id == stuck_task_id)
    ).scalar_one_or_none()
    celery_meta_before = db_session_sync.execute(
        select(CeleryTaskMeta).where(CeleryTaskMeta.task_id == stuck_task_id)
    ).scalar_one_or_none()
    assert task_history_before is not None
    assert celery_meta_before is not None

    result = cleanup_stuck_tasks_task(
        stuck_pending_hours=stuck_pending_hours,
        stuck_started_hours=stuck_started_hours,
    )

    assert result["status"] == "completed"
    assert result["number_of_stuck_pending_deleted"] == 1

    task_history_after = db_session_sync.execute(
        select(TaskHistoryDB).where(TaskHistoryDB.task_id == stuck_task_id)
    ).scalar_one_or_none()
    celery_meta_after = db_session_sync.execute(
        select(CeleryTaskMeta).where(CeleryTaskMeta.task_id == stuck_task_id)
    ).scalar_one_or_none()
    assert task_history_after is None
    assert celery_meta_after is None
