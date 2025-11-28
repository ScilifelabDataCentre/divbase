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


async def get_tasks_for_user_pg(
    db: AsyncSession,
    user_id: int,
    is_admin: bool = False,
) -> list[dict]:
    """
    Fetch task data (CeleryTaskMeta + TaskHistoryDB) for a user in a single query.
    Filters on user_id.
    """
    stmt = select(
        *CeleryTaskMeta.__table__.c,
        TaskHistoryDB.created_at,
        TaskHistoryDB.started_at,
        TaskHistoryDB.completed_at,
    ).join(TaskHistoryDB, CeleryTaskMeta.task_id == TaskHistoryDB.task_id)

    if not is_admin:
        stmt = stmt.where(TaskHistoryDB.user_id == user_id)

    result = await db.execute(stmt)
    rows = result.fetchall()
    return [dict(row._mapping) for row in rows]


async def get_tasks_for_user_and_project_pg(
    db: AsyncSession,
    user_id: int,
    project_id: int,
    is_admin: bool = False,
) -> list[dict]:
    """
    Fetch task data (CeleryTaskMeta + TaskHistoryDB) for a user and project in a single query.
    Filters on user_id and project_id.
    """
    stmt = (
        select(
            *CeleryTaskMeta.__table__.c,
            TaskHistoryDB.created_at,
            TaskHistoryDB.started_at,
            TaskHistoryDB.completed_at,
        )
        .join(TaskHistoryDB, CeleryTaskMeta.task_id == TaskHistoryDB.task_id)
        .where(TaskHistoryDB.project_id == project_id)
    )

    if not is_admin:
        stmt = stmt.where(TaskHistoryDB.user_id == user_id)

    result = await db.execute(stmt)
    rows = result.fetchall()
    return [dict(row._mapping) for row in rows]


async def get_tasks_for_project_pg(
    db: AsyncSession,
    project_id: int,
) -> list[dict]:
    """
    Fetch task data (CeleryTaskMeta + TaskHistoryDB) for a project in a single query.
    API layer handles permission (MANAGE user role required).
    """
    stmt = (
        select(
            *CeleryTaskMeta.__table__.c,
            TaskHistoryDB.created_at,
            TaskHistoryDB.started_at,
            TaskHistoryDB.completed_at,
        )
        .join(TaskHistoryDB, CeleryTaskMeta.task_id == TaskHistoryDB.task_id)
        .where(TaskHistoryDB.project_id == project_id)
    )

    result = await db.execute(stmt)
    rows = result.fetchall()
    return [dict(row._mapping) for row in rows]


async def get_task_by_id_if_user_allowed(
    db: AsyncSession,
    task_id: str,
    user_id: int,
    is_admin: bool,
) -> dict | None:
    """
    Fetch task by ID only if user has permission to view it.
    Returns task dict if found and allowed, None otherwise.

    This combines permission check and data fetch in a single query.
    """
    # TODO move the permissions check up to the routes layer?
    stmt = (
        select(
            *CeleryTaskMeta.__table__.c,
            TaskHistoryDB.created_at,
            TaskHistoryDB.started_at,
            TaskHistoryDB.completed_at,
        )
        .join(TaskHistoryDB, CeleryTaskMeta.task_id == TaskHistoryDB.task_id)
        .where(TaskHistoryDB.task_id == task_id)
    )

    if is_admin:
        pass
    else:
        # When not admin, must be submitter of Manager user to view task
        stmt = stmt.where(
            (TaskHistoryDB.user_id == user_id)
            | (
                TaskHistoryDB.project_id.in_(
                    select(ProjectMembershipDB.project_id).where(
                        ProjectMembershipDB.user_id == user_id, ProjectMembershipDB.role == ProjectRoles.MANAGE
                    )
                )
            )
        )

    result = await db.execute(stmt)
    row = result.first()

    return dict(row._mapping) if row else None
