"""
CRUD operations for task history.
"""

import logging

from sqlalchemy import join, select
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.models.projects import ProjectDB, ProjectMembershipDB, ProjectRoles
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


async def get_allowed_task_ids_for_user(
    db: AsyncSession,
    user_id: int,
    is_admin: bool,
) -> set[str]:
    """
    Get all task IDs the user has permission to view.

    Task IDs are unique keys in the TaskHistoryDB table, so a set is not strictly needed. But a set is used for faster lookup.
    """
    if is_admin:
        stmt = select(TaskHistoryDB.task_id)
    else:
        stmt = select(TaskHistoryDB.task_id).where(TaskHistoryDB.user_id == user_id)
    result = await db.execute(stmt)
    rows = result.fetchall()
    allowed_task_ids = set()
    for row in rows:
        allowed_task_ids.add(row[0])
    return allowed_task_ids


async def get_allowed_task_ids_for_user_and_project(
    db: AsyncSession,
    user_id: int,
    project_id: int,
    is_admin: bool,
) -> set[str]:
    """
    Get all task IDs the user has permission to view.

    Task IDs are unique keys in the TaskHistoryDB table, so a set is not strictly needed. But a set is used for faster lookup.
    """
    if is_admin:
        stmt = select(TaskHistoryDB.task_id).where(TaskHistoryDB.project_id == project_id)
    else:
        stmt = (
            select(TaskHistoryDB.task_id)
            .where(TaskHistoryDB.user_id == user_id)
            .where(TaskHistoryDB.project_id == project_id)
        )
    result = await db.execute(stmt)
    rows = result.fetchall()
    allowed_task_ids = set()
    for row in rows:
        allowed_task_ids.add(row[0])
    return allowed_task_ids


async def filter_task_ids_by_project_name(db: AsyncSession, task_ids: set[str], project_name: str) -> set[str]:
    """
    Filter a Set of task IDs by project name. Uses relationship between TaskHistoryDB and ProjectDB to make a single query with a table join.

    Task IDs are unique keys in the TaskHistoryDB table, so a set is not strictly needed. But a set is used for faster lookup.
    """

    stmt = (
        select(TaskHistoryDB.task_id)
        .select_from(join(TaskHistoryDB, ProjectDB, TaskHistoryDB.project_id == ProjectDB.id))
        .where(ProjectDB.name == project_name, TaskHistoryDB.task_id.in_(task_ids))
    )
    result = await db.execute(stmt)

    rows = result.fetchall()
    allowed_task_ids = set()
    for row in rows:
        allowed_task_ids.add(row[0])
    return allowed_task_ids


async def get_allowed_task_ids_for_project(
    db: AsyncSession,
    project_id: int,
) -> set[str]:
    """
    Get all task IDs for a project.

    Task IDs are unique keys in the TaskHistoryDB table, so a set is not strictly needed. But a set is used for faster lookup.
    """

    stmt = select(TaskHistoryDB.task_id).where(TaskHistoryDB.project_id == project_id)
    result = await db.execute(stmt)
    rows = result.fetchall()
    allowed_task_ids = set()
    for row in rows:
        allowed_task_ids.add(row[0])
    return allowed_task_ids


async def check_user_can_view_task_id(
    db: AsyncSession,
    task_id: str,
    user_id: int,
    is_admin: bool,
) -> bool:
    """
    Check if a user has permission to view the task history for a specific task ID.

    Admin is allowed to view all task IDs.
    Manage user is allowed to see all task IDs for the project the they manage.
    Edit user is only allowed see task IDs for tasks they submitted themselves.

    Uses explicit joins the two tables to be joined share a foreign key to a third table,
    but no direct FK between them.
    """

    if is_admin:
        return True

    # Case: submitting user
    stmt = select(TaskHistoryDB.task_id).where(
        TaskHistoryDB.task_id == task_id,
        TaskHistoryDB.user_id == user_id,
    )
    result = await db.execute(stmt)
    if result.scalar_one_or_none():
        return True

    # Case: manager user
    stmt = (
        select(TaskHistoryDB.task_id)
        .join(ProjectMembershipDB, TaskHistoryDB.project_id == ProjectMembershipDB.project_id)
        .where(
            TaskHistoryDB.task_id == task_id,
            ProjectMembershipDB.user_id == user_id,
            ProjectMembershipDB.role == ProjectRoles.MANAGE,
        )
    )
    result = await db.execute(stmt)
    return result.scalar_one_or_none() is not None
