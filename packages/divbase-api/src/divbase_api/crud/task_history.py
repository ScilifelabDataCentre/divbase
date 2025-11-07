"""
CRUD operations for task history.
"""

import logging

from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.models.task_history import TaskHistoryDB, TaskStatus

logger = logging.getLogger(__name__)


async def record_pending_task(db: AsyncSession, task_id: str, user_id: int, project_id: int):
    """
    Record a new pending task in the TaskHistoryDB table. Intended to be run when task is submitted to queue.
    """
    entry = TaskHistoryDB(
        task_id=str(task_id),
        user_id=user_id,
        project_id=project_id,
        status=TaskStatus.PENDING,
    )
    db.add(entry)
    await db.commit()
    await db.refresh(entry)


async def check_user_can_view_task_history(
    db: AsyncSession,
    task_id: str,
    user_id: int,
    is_admin: bool,
) -> bool:
    """
    Check if the user has permission to view this task.
    Returns False if task not found in postgres TaskHistoryDB table or in results backend
    (redis backend may purge old tasks).
    """
    if is_admin:
        return True

    try:
        stmt = select(TaskHistoryDB).where(TaskHistoryDB.task_id == task_id)
        result = await db.execute(stmt)
        entry = result.scalar_one_or_none()
        if entry:
            return entry.user_id == user_id
        else:
            return False
    except Exception as e:
        logger.error(f"Database error checking permission for task {task_id}: {e}")
        return False
