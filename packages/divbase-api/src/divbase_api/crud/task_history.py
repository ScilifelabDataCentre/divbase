"""
CRUD operations for task history.
"""

import logging

from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.models.projects import ProjectMembershipDB, ProjectRoles
from divbase_api.models.task_history import CeleryTaskMeta, TaskHistoryDB, TaskStartedAtDB
from divbase_api.models.users import UserDB

logger = logging.getLogger(__name__)


async def get_tasks_pg(
    db: AsyncSession,
    user_id: int | None = None,
    project_id: int | None = None,
    user_task_id: int | None = None,
    is_admin: bool = False,
    require_manager_role: bool = False,
) -> list[dict] | None:
    """
    Dynamic crud to fetch task history (CeleryTaskMeta + TaskHistoryDB) with different filters.
    If task_id is provided, returns a single dict or None.
    Otherwise, returns a list of dicts.

    System tasks (project_id=NULL) are excluded by default.

    Note: Packing of results into a pydantic model happens in during deserialization in the service layer.
    """
    stmt = (
        select(
            *CeleryTaskMeta.__table__.c,
            TaskHistoryDB.id.label("user_task_id"),
            TaskHistoryDB.created_at,
            TaskStartedAtDB.started_at,
            UserDB.email.label(
                "submitter_email"
            ),  # TODO consider not using label here since email is unique in this selection
        )
        .join(TaskHistoryDB, CeleryTaskMeta.task_id == TaskHistoryDB.task_id)
        .join(UserDB, TaskHistoryDB.user_id == UserDB.id)
        .join(TaskStartedAtDB, TaskHistoryDB.task_id == TaskStartedAtDB.task_id)
    )

    if user_task_id:
        stmt = stmt.where(TaskHistoryDB.id == user_task_id)
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
    if user_task_id:
        row = result.first()
        return [dict(row._mapping)] if row else None
    else:
        rows = result.fetchall()
        return [dict(row._mapping) for row in rows]


async def create_task_history_entry(
    db: AsyncSession,
    task_id: str,
    user_id: int,
    project_id: int,
) -> int:
    """
    Create a new task history entry. Returns the primary key ID from the table, which is the user's job id (but different from celery task_id).
    """
    task_history_entry = TaskHistoryDB(task_id=task_id, user_id=user_id, project_id=project_id)
    db.add(task_history_entry)
    await db.commit()
    await db.refresh(task_history_entry)
    return task_history_entry.id
