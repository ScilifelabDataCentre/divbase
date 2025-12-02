"""
CRUD operations for task history.
"""

import logging

from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.models.projects import ProjectMembershipDB, ProjectRoles
from divbase_api.models.task_history import CeleryTaskMeta, TaskHistoryDB, TaskStatus

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


async def get_tasks_pg(
    db: AsyncSession,
    user_id: int | None = None,
    project_id: int | None = None,
    task_id: str | None = None,
    is_admin: bool = False,
    require_manager_role: bool = False,
) -> list[dict] | dict | None:
    """
    Dynamic crud to fetch task history (CeleryTaskMeta + TaskHistoryDB) with different filters.
    If task_id is provided, returns a single dict or None.
    Otherwise, returns a list of dicts.

    Note: Packing of results into a pydantic model happens in during deserialization in the service layer.
    """
    stmt = select(
        *CeleryTaskMeta.__table__.c,
        TaskHistoryDB.created_at,
        TaskHistoryDB.started_at,
        TaskHistoryDB.completed_at,
    ).join(TaskHistoryDB, CeleryTaskMeta.task_id == TaskHistoryDB.task_id)

    if task_id:
        stmt = stmt.where(TaskHistoryDB.task_id == task_id)
    if project_id is not None:
        stmt = stmt.where(TaskHistoryDB.project_id == project_id)
    if user_id is not None and not is_admin:
        if require_manager_role:
            stmt = stmt.where(
                (TaskHistoryDB.user_id == user_id)
                | (
                    TaskHistoryDB.project_id.in_(
                        select(ProjectMembershipDB.project_id).where(
                            ProjectMembershipDB.user_id == user_id,
                            ProjectMembershipDB.role == ProjectRoles.MANAGE,
                        )
                    )
                )
            )
        else:
            stmt = stmt.where(TaskHistoryDB.user_id == user_id)

    result = await db.execute(stmt)
    if task_id:
        row = result.first()
        return dict(row._mapping) if row else None
    else:
        rows = result.fetchall()
        return [dict(row._mapping) for row in rows]
