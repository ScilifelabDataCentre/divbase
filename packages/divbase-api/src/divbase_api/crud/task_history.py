"""
CRUD operations for task history..
"""

import logging

from divbase_api.models.task_history import TaskHistoryDB, TaskStatus
from divbase_api.worker.worker_db import SyncSessionLocal

logger = logging.getLogger(__name__)


def record_pending_task(task_id: str, user_id: int, project_id: int):
    with SyncSessionLocal() as session:
        metadata = TaskHistoryDB(
            task_id=str(task_id), user_id=user_id, project_id=project_id, status=TaskStatus.PENDING
        )
        session.add(metadata)
        session.commit()
