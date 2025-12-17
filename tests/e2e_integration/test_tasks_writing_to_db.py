"""
Test for the postgres task history upsert functionality in _upsert_task_timestamps.
Simulates multipls threads of workers using threading library.
"""

# import threading
# from unittest.mock import patch

# import pytest
# from sqlalchemy import select, text

# from divbase_api.models.task_history import TaskHistoryDB
# from divbase_api.services.s3_client import create_s3_file_manager
# from divbase_api.worker.tasks import _upsert_task_timestamps, update_vcf_dimensions_task
# from divbase_api.worker.worker_db import SyncSessionLocal
# from divbase_api.models import TaskStartedAtDB

# # @pytest.fixture
# def clean_task_history(db_session_sync):
#     """Remove all task history entries before each test."""
#     db_session_sync.query(TaskHistoryDB).delete()
#     db_session_sync.commit()
#     yield
#     db_session_sync.query(TaskHistoryDB).delete()
#     db_session_sync.commit()


# @pytest.fixture(autouse=True)
# def clean_celery_taskmeta(db_session_sync):
#     """Remove all celery_taskmeta entries before each test."""
#     db_session_sync.execute(text("DELETE FROM celery_taskmeta"))
#     db_session_sync.commit()
#     yield
#     db_session_sync.execute(text("DELETE FROM celery_taskmeta"))
#     db_session_sync.commit()


# def test_insert_new_beat_task(clean_task_history, clean_task_started_at):
#     """Test inserting a new beat-scheduled task (no existing entry in db prior to execution)."""
#     task_id = "beat-task-123"

#     # Simulate the after_task_publish signal creating TaskHistoryDB entry
#     with SyncSessionLocal() as db:
#         entry = TaskHistoryDB(task_id=task_id, user_id=None, project_id=None)
#         db.add(entry)
#         db.commit()

#     # Simulate the task_prerun signal creating TaskStartedAtDB entry
#     with SyncSessionLocal() as db:
#         started_entry = TaskStartedAtDB(task_id=task_id)
#         db.add(started_entry)
#         db.commit()

#     # Verify TaskHistoryDB entry
#     with SyncSessionLocal() as db:
#         entry = db.execute(select(TaskHistoryDB).where(TaskHistoryDB.task_id == task_id)).scalar_one()
#         assert entry.task_id == task_id
#         assert entry.user_id is None
#         assert entry.project_id is None

#     # Verify TaskStartedAtDB entry
#     with SyncSessionLocal() as db:
#         started = db.execute(select(TaskStartedAtDB).where(TaskStartedAtDB.task_id == task_id)).scalar_one()
#         assert started.task_id == task_id
#         assert started.started_at is not None


# def test_update_existing_api_task(clean_task_history, clean_task_started_at):
#     """Test that task_prerun signal adds started_at for an existing API-submitted task."""
#     task_id = "api-task-456"
#     user_id = 1
#     project_id = 1

#     with SyncSessionLocal() as db:
#         entry = TaskHistoryDB(
#             task_id=task_id,
#             user_id=user_id,
#             project_id=project_id,
#         )
#         db.add(entry)
#         db.commit()

#     # Simulate task_prerun signal creating TaskStartedAtDB entry
#     with SyncSessionLocal() as db:
#         started_entry = TaskStartedAtDB(task_id=task_id)
#         db.add(started_entry)
#         db.commit()

#     # Verify TaskHistoryDB entry is unchanged
#     with SyncSessionLocal() as db:
#         entry = db.execute(select(TaskHistoryDB).where(TaskHistoryDB.task_id == task_id)).scalar_one()
#         assert entry.user_id == user_id
#         assert entry.project_id == project_id

#     # Verify TaskStartedAtDB entry exists
#     with SyncSessionLocal() as db:
#         started = db.execute(select(TaskStartedAtDB).where(TaskStartedAtDB.task_id == task_id)).scalar_one()
#         assert started.started_at is not None


# def test_concurrent_inserts_race_condition(clean_task_history):
#     """
#     Test that concurrent inserts don't create duplicates.
#     Simulates two workers trying to create the same entry simultaneously in TaskHistoryDB.
#     Could happen if two Celery signals fire at the same time for the same task, e.g due to retries, race conditions, or worker thread pools.
#     """
#     task_id = "race-task-789"

#     errors = []

#     def worker_insert():
#         try:
#             _upsert_task_timestamps(
#                 task_id=task_id,
#                 set_started_at=True,
#             )
#         except Exception as e:
#             errors.append(e)

#     # Start two threads simultaneously
#     thread1 = threading.Thread(target=worker_insert)
#     thread2 = threading.Thread(target=worker_insert)

#     thread1.start()
#     thread2.start()

#     thread1.join()
#     thread2.join()

#     assert len(errors) == 0

#     with SyncSessionLocal() as db:
#         entries = db.execute(select(TaskHistoryDB).where(TaskHistoryDB.task_id == task_id)).all()
#         assert len(entries) == 1


# def test_idempotent_updates(clean_task_history):
#     """Test that multiple updates are idempotent, in this case that started_at should still be the original value after upsert."""
#     task_id = "idempotent-task-111"

#     _upsert_task_timestamps(
#         task_id=task_id,
#         set_started_at=True,
#     )

#     with SyncSessionLocal() as db:
#         first_started_at = db.execute(
#             select(TaskHistoryDB.started_at).where(TaskHistoryDB.task_id == task_id)
#         ).scalar_one()

#     _upsert_task_timestamps(
#         task_id=task_id,
#         set_completed_at=True,
#     )

#     with SyncSessionLocal() as db:
#         entry = db.execute(select(TaskHistoryDB).where(TaskHistoryDB.task_id == task_id)).scalar_one()

#         assert entry.started_at == first_started_at
#         assert entry.completed_at is not None


# def test_no_duplicate_entries_created(clean_task_history):
#     """Verify no duplicate task_id entries can exist in TaskHistoryDB."""
#     task_id = "duplicate-check-222"

#     _upsert_task_timestamps(
#         task_id=task_id,
#         set_started_at=True,
#     )

#     _upsert_task_timestamps(
#         task_id=task_id,
#         set_completed_at=True,
#     )

#     with SyncSessionLocal() as db:
#         count = db.execute(select(TaskHistoryDB).where(TaskHistoryDB.task_id == task_id)).all()
#         assert len(count) == 1


# def test_database_connection_error_handling(clean_task_history):
#     """Test that database errors are raised."""
#     task_id = "error-task-333"

#     with patch("divbase_api.worker.tasks.SyncSessionLocal") as mock_session:
#         mock_session.return_value.__enter__.return_value.execute.side_effect = Exception("DB Error")

#         with pytest.raises(Exception, match="DB Error"):
#             _upsert_task_timestamps(
#                 task_id=task_id,
#             )


# def test_concurrent_inserts_race_condition_real_task(clean_task_history, CONSTANTS, project_map):
#     """
#     Test concurrent execution of real Celery tasks to verify no duplicate task history entries.
#     Submit same task multiple times concurrently to verify that there are no duplicate task history entries.
#     """

#     project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
#     bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
#     project_id = project_map[project_name]

#     task_ids = []
#     async_results = []

#     with patch("divbase_api.worker.tasks.create_s3_file_manager") as mock_create_s3_manager:
#         mock_create_s3_manager.side_effect = lambda url=None: create_s3_file_manager(url=CONSTANTS["MINIO_URL"])

#         for _ in range(3):
#             result = update_vcf_dimensions_task.apply_async(
#                 kwargs={
#                     "bucket_name": bucket_name,
#                     "project_id": project_id,
#                     "project_name": project_name,
#                     "user_id": 1,
#                 }
#             )
#             task_ids.append(str(result.id))
#             async_results.append(result)

#         for result in async_results:
#             result.get(timeout=30)

#     # Verify each task_id has exactly one entry in TaskHistoryDB
#     with SyncSessionLocal() as db:
#         for tid in task_ids:
#             entries = db.execute(select(TaskHistoryDB).where(TaskHistoryDB.task_id == tid)).all()
#             assert len(entries) == 1, f"Task {tid} has {len(entries)} entries, expected 1"
