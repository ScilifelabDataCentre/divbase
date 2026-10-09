"""
Test that the task_prerun signal handler flags a task that was started before (i.e. has been redelivered).
"""

import uuid
from types import SimpleNamespace

from sqlalchemy import delete

from divbase_api.models.task_history import TaskStartedAtDB
from divbase_api.worker.tasks import handle_task_started


def test_prerun_handler_flags_task_started_before(db_session_sync):
    task_id = f"redelivery-test-{uuid.uuid4()}"
    first_run = SimpleNamespace(request=SimpleNamespace())
    redelivered_run = SimpleNamespace(request=SimpleNamespace())

    try:
        handle_task_started(sender=first_run, task_id=task_id)
        handle_task_started(sender=redelivered_run, task_id=task_id)

        assert first_run.request.task_already_ran is False
        assert redelivered_run.request.task_already_ran is True
    finally:
        db_session_sync.execute(delete(TaskStartedAtDB).where(TaskStartedAtDB.task_id == task_id))
        db_session_sync.commit()


def test_prerun_handler_does_not_flag_different_tasks(db_session_sync):
    task_ids = [f"redelivery-test-{uuid.uuid4()}" for _ in range(2)]
    runs = [SimpleNamespace(request=SimpleNamespace()) for _ in task_ids]

    try:
        for run, task_id in zip(runs, task_ids, strict=True):
            handle_task_started(sender=run, task_id=task_id)

        assert all(run.request.task_already_ran is False for run in runs)
    finally:
        db_session_sync.execute(delete(TaskStartedAtDB).where(TaskStartedAtDB.task_id.in_(task_ids)))
        db_session_sync.commit()
