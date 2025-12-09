"""
Test for the postgres task history upsert functionality in _update_task_status_in_pg.
Simulates multipls threads of workers using threading library.

Several of these tests are unit test that need to be in e2e_integration folder since it needs postgres from the testing stack to test against.
"""

import threading
from unittest.mock import patch

import pytest
from sqlalchemy import select, text

from divbase_api.models.task_history import TaskHistoryDB, TaskStatus
from divbase_api.services.s3_client import create_s3_file_manager
from divbase_api.worker.tasks import _update_task_status_in_pg, task_pending_handler, update_vcf_dimensions_task
from divbase_api.worker.worker_db import SyncSessionLocal


@pytest.fixture
def clean_task_history(db_session_sync):
    """
    Remove all task history entries before each test. Use raw SQL for celery_taskmeta table since it is created and managed by Celery and not by the SQLAlchemy ORM.
    """
    db_session_sync.query(TaskHistoryDB).delete()
    db_session_sync.commit()
    db_session_sync.execute(text("DELETE FROM celery_taskmeta"))
    db_session_sync.commit()


def test_insert_new_beat_task(clean_task_history):
    """Test inserting a new beat-scheduled task (no existing entry in db prior to execution of update function)."""
    task_id = "beat-task-123"

    _update_task_status_in_pg(
        task_id=task_id,
        status=TaskStatus.STARTED,
        set_started_at=True,
    )

    with SyncSessionLocal() as db:
        entry = db.execute(select(TaskHistoryDB).where(TaskHistoryDB.task_id == task_id)).scalar_one()

        assert entry.task_id == task_id
        assert entry.status == TaskStatus.STARTED
        assert entry.user_id is None
        assert entry.project_id is None
        assert entry.started_at is not None
        assert entry.completed_at is None


def test_update_existing_api_task(clean_task_history):
    """Test updating an existing API-submitted task (entry exists in db prior to execution of update function)."""
    task_id = "api-task-456"
    user_id = 1
    project_id = 1

    with SyncSessionLocal() as db:
        entry = TaskHistoryDB(
            task_id=task_id,
            user_id=user_id,
            project_id=project_id,
            status=TaskStatus.PENDING,
        )
        db.add(entry)
        db.commit()

    _update_task_status_in_pg(
        task_id=task_id,
        status=TaskStatus.STARTED,
        set_started_at=True,
    )

    with SyncSessionLocal() as db:
        entry = db.execute(select(TaskHistoryDB).where(TaskHistoryDB.task_id == task_id)).scalar_one()
        assert entry.user_id == user_id
        assert entry.project_id == project_id
        assert entry.status == TaskStatus.STARTED
        assert entry.started_at is not None


def test_concurrent_inserts_race_condition(clean_task_history):
    """
    Test that concurrent inserts don't create duplicates.
    Simulates two workers trying to create the same entry simultaneously in TaskHistoryDB.
    Could happen if two Celery signals fire at the same time for the same task, e.g due to retries, race conditions, or worker thread pools.
    """
    task_id = "race-task-789"

    errors = []

    def worker_insert():
        try:
            _update_task_status_in_pg(
                task_id=task_id,
                status=TaskStatus.STARTED,
                set_started_at=True,
            )
        except Exception as e:
            errors.append(e)

    # Start two threads simultaneously
    thread1 = threading.Thread(target=worker_insert)
    thread2 = threading.Thread(target=worker_insert)

    thread1.start()
    thread2.start()

    thread1.join()
    thread2.join()

    assert len(errors) == 0

    with SyncSessionLocal() as db:
        entries = db.execute(select(TaskHistoryDB).where(TaskHistoryDB.task_id == task_id)).all()
        assert len(entries) == 1


def test_update_with_error_message(clean_task_history):
    """Test updating task with error message."""
    task_id = "error-task-999"
    error_msg = "Something went wrong"

    _update_task_status_in_pg(
        task_id=task_id,
        status=TaskStatus.STARTED,
        set_started_at=True,
    )

    _update_task_status_in_pg(
        task_id=task_id,
        status=TaskStatus.FAILURE,
        error_msg=error_msg,
        set_completed_at=True,
    )

    with SyncSessionLocal() as db:
        entry = db.execute(select(TaskHistoryDB).where(TaskHistoryDB.task_id == task_id)).scalar_one()

        assert entry.status == TaskStatus.FAILURE
        assert entry.error_message == error_msg
        assert entry.completed_at is not None


def test_idempotent_status_updates(clean_task_history):
    """Test that multiple updates are idempotent, in this case that started_at should still be the original value after upsert."""
    task_id = "idempotent-task-111"

    _update_task_status_in_pg(
        task_id=task_id,
        status=TaskStatus.STARTED,
        set_started_at=True,
    )

    with SyncSessionLocal() as db:
        first_started_at = db.execute(
            select(TaskHistoryDB.started_at).where(TaskHistoryDB.task_id == task_id)
        ).scalar_one()

    _update_task_status_in_pg(
        task_id=task_id,
        status=TaskStatus.SUCCESS,
        set_completed_at=True,
    )

    with SyncSessionLocal() as db:
        entry = db.execute(select(TaskHistoryDB).where(TaskHistoryDB.task_id == task_id)).scalar_one()

        assert entry.started_at == first_started_at
        assert entry.status == TaskStatus.SUCCESS
        assert entry.completed_at is not None


def test_no_duplicate_entries_created(clean_task_history):
    """Verify no duplicate task_id entries can exist in TaskHistoryDB."""
    task_id = "duplicate-check-222"

    _update_task_status_in_pg(
        task_id=task_id,
        status=TaskStatus.STARTED,
        set_started_at=True,
    )

    _update_task_status_in_pg(
        task_id=task_id,
        status=TaskStatus.SUCCESS,
        set_completed_at=True,
    )

    with SyncSessionLocal() as db:
        count = db.execute(select(TaskHistoryDB).where(TaskHistoryDB.task_id == task_id)).all()
        assert len(count) == 1


def test_database_connection_error_handling(clean_task_history):
    """Test that database errors are logged but don't crash."""
    task_id = "error-task-333"

    with (
        patch("divbase_api.worker.tasks.SyncSessionLocal") as mock_session,
        patch("divbase_api.worker.tasks.logger") as mock_logger,
    ):
        mock_session.return_value.__enter__.return_value.execute.side_effect = Exception("DB Error")

        # Should not raise
        _update_task_status_in_pg(
            task_id=task_id,
            status=TaskStatus.STARTED,
        )

        # Verify error was logged
        mock_logger.error.assert_called_once()
        assert "Failed to update task status" in str(mock_logger.error.call_args)


def test_concurrent_inserts_race_condition_real_task(clean_task_history, CONSTANTS, project_map):
    """
    Test concurrent execution of real Celery tasks to verify no duplicate task history entries.
    Submit same task multiple times concurrently to verify that there are no duplicate task history entries.
    """

    project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    project_id = project_map[project_name]
    user_id = 1

    task_ids = []
    async_results = []

    with patch("divbase_api.worker.tasks.create_s3_file_manager") as mock_create_s3_manager:
        mock_create_s3_manager.side_effect = lambda url=None: create_s3_file_manager(url=CONSTANTS["MINIO_URL"])

        for _ in range(3):
            result = update_vcf_dimensions_task.apply_async(
                kwargs={
                    "bucket_name": bucket_name,
                    "project_id": project_id,
                    "project_name": project_name,
                    "user_id": user_id,
                }
            )
            task_ids.append(str(result.id))
            async_results.append(result)

        for result in async_results:
            result.get(timeout=30)

    # Verify each task_id has exactly one entry in TaskHistoryDB
    with SyncSessionLocal() as db:
        for tid in task_ids:
            entries = db.execute(select(TaskHistoryDB).where(TaskHistoryDB.task_id == tid)).all()
            assert len(entries) == 1, f"Task {tid} has {len(entries)} entries, expected 1"


def test_pending_status_created_on_task_publish(clean_task_history):
    """Test that PENDING status is created when task is published (after_task_publish signal)."""
    task_id = "pending-task-111"
    user_id = 1
    project_id = 2

    # Simulate what happens in after_task_publish signal
    headers = {"id": task_id}
    body = (
        [],  # args
        {"user_id": user_id, "project_id": project_id, "bucket_name": "test-bucket"},  # kwargs
        {"callbacks": None, "errbacks": None, "chain": None, "chord": None},  # embed
    )

    task_pending_handler(sender=None, headers=headers, body=body)

    with SyncSessionLocal() as db:
        entry = db.execute(select(TaskHistoryDB).where(TaskHistoryDB.task_id == task_id)).scalar_one()
        assert entry.task_id == task_id
        assert entry.status == TaskStatus.PENDING
        assert entry.user_id == user_id
        assert entry.project_id == project_id
        assert entry.started_at is None
        assert entry.completed_at is None


def test_pending_status_without_user_project(clean_task_history):
    """Test that PENDING status works for tasks without user_id/project_id (=beat tasks)."""
    task_id = "beat-pending-task-222"

    headers = {"id": task_id}
    body = (
        [],  # args
        {},  # kwargs without user_id/project_id
        {},  # embed
    )

    task_pending_handler(sender=None, headers=headers, body=body)

    with SyncSessionLocal() as db:
        entry = db.execute(select(TaskHistoryDB).where(TaskHistoryDB.task_id == task_id)).scalar_one()

        assert entry.status == TaskStatus.PENDING
        assert entry.user_id is None
        assert entry.project_id is None


def test_pending_then_started_preserves_user_project(clean_task_history):
    """
    Test that PENDING -> STARTED transition preserves user_id and project_id.
    Checks that the ON CONFLICT... UPDATE logic works as intended.
    """
    task_id = "transition-task-333"
    user_id = 1
    project_id = 2

    headers = {"id": task_id}
    body = ([], {"user_id": user_id, "project_id": project_id}, {})

    task_pending_handler(sender=None, headers=headers, body=body)

    _update_task_status_in_pg(
        task_id=task_id,
        status=TaskStatus.STARTED,
        set_started_at=True,
    )

    with SyncSessionLocal() as db:
        entry = db.execute(select(TaskHistoryDB).where(TaskHistoryDB.task_id == task_id)).scalar_one()

        assert entry.status == TaskStatus.STARTED
        assert entry.user_id == user_id
        assert entry.project_id == project_id
        assert entry.started_at is not None


def test_pending_status_concurrent_publish(clean_task_history):
    """
    Test that concurrent after_task_publish calls don't create duplicates.
    Checks that the idempotency of the UPSERT works as intended.
    """
    task_id = "concurrent-pending-666"
    user_id = 1
    project_id = 2

    headers = {"id": task_id}
    body = ([], {"user_id": user_id, "project_id": project_id}, {})

    errors = []

    def publish_task():
        try:
            task_pending_handler(sender=None, headers=headers, body=body)
        except Exception as e:
            errors.append(e)

    # Simulate concurrent publishes
    thread1 = threading.Thread(target=publish_task)
    thread2 = threading.Thread(target=publish_task)

    thread1.start()
    thread2.start()

    thread1.join()
    thread2.join()

    assert len(errors) == 0

    # Verify only one entry exists
    with SyncSessionLocal() as db:
        entries = db.execute(select(TaskHistoryDB).where(TaskHistoryDB.task_id == task_id)).all()
        assert len(entries) == 1


def test_full_lifecycle_with_pending(clean_task_history):
    """
    Test complete task lifecycle: PENDING -> STARTED -> SUCCESS.
    This can be trick to test with a real task since the step from PENDING -> STARTED will be instataneous when there are idle workers.
    """
    task_id = "lifecycle-task-777"
    user_id = 1
    project_id = 2

    # Step 1: PENDING (from after_task_publish)
    headers = {"id": task_id}
    body = ([], {"user_id": user_id, "project_id": project_id}, {})

    task_pending_handler(sender=None, headers=headers, body=body)

    with SyncSessionLocal() as db:
        entry = db.execute(select(TaskHistoryDB).where(TaskHistoryDB.task_id == task_id)).scalar_one()
        assert entry.status == TaskStatus.PENDING
        assert entry.started_at is None
        assert entry.completed_at is None

    # Step 2: STARTED (from task_prerun)
    _update_task_status_in_pg(task_id=task_id, status=TaskStatus.STARTED, set_started_at=True)

    with SyncSessionLocal() as db:
        entry = db.execute(select(TaskHistoryDB).where(TaskHistoryDB.task_id == task_id)).scalar_one()
        assert entry.status == TaskStatus.STARTED
        assert entry.started_at is not None
        assert entry.completed_at is None

    # Step 3: SUCCESS (from task_success)
    _update_task_status_in_pg(task_id=task_id, status=TaskStatus.SUCCESS, set_completed_at=True)

    with SyncSessionLocal() as db:
        entry = db.execute(select(TaskHistoryDB).where(TaskHistoryDB.task_id == task_id)).scalar_one()
        assert entry.status == TaskStatus.SUCCESS
        assert entry.user_id == user_id
        assert entry.project_id == project_id
        assert entry.started_at is not None
        assert entry.completed_at is not None
