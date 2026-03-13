"""
CRUD operations for task history.
"""

import logging
from datetime import datetime, timedelta, timezone

from sqlalchemy import and_, delete, or_, select, update
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.models.projects import ProjectMembershipDB, ProjectRoles
from divbase_api.models.task_history import CeleryTaskMeta, TaskHistoryDB, TaskStartedAtDB
from divbase_api.models.users import UserDB
from divbase_api.worker.tasks import TaskName

logger = logging.getLogger(__name__)

ACTIVE_CELERY_STATUSES = {"PENDING", "STARTED", "RETRY"}


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
    user_id: int,
    project_id: int,
    task_id: str | None = None,
) -> int:
    """
    Create a new task history entry. Returns the primary key ID from the table, which is the DivBase job ID
    (which is an integer rather than the celery task_id which is a UUID).

    The task_id parameter is optional to allow creating entries before the celery task is created (e.g., for bcftools pipe tasks).

    """
    task_history_entry = TaskHistoryDB(task_id=task_id, user_id=user_id, project_id=project_id)
    db.add(task_history_entry)
    await db.commit()
    await db.refresh(task_history_entry)
    return task_history_entry.id


async def update_task_history_entry_with_celery_task_id(
    db: AsyncSession,
    job_id: int,
    task_id: str | None = None,
) -> None:
    """
    Updates an existing task history entry with the celery task ID. The intended use case is to first create a table entry without a celery task ID using
    create_task_history_entry(), then submit a celery task with .apply_async and get the celery task ID, and then update the table entry with the celery task ID using this function.
    """

    stmt = (
        update(TaskHistoryDB)
        .where(
            TaskHistoryDB.id == job_id,
            TaskHistoryDB.task_id.is_(None),
        )
        .values(task_id=task_id)
    )
    result = await db.execute(stmt)
    if result.rowcount == 0:
        raise ValueError(f"TaskHistoryDB entry with id={job_id} not found or already finalized with task_id")
    await db.commit()


async def create_provisional_dimensions_reservation(db: AsyncSession, user_id: int, project_id: int) -> int:
    """
    Creates a provisional reservation for dimensions update by creating a TaskHistoryDB entry with task_id=None and
    task_name="update_vcf_dimensions_task". This allows the update_vcf_dimensions_task to be associated with this entry
    once the celery task is created and gets its task_id, and allows the task history to be visible in the UI as soon
    as the reservation is made rather than waiting for the celery task to be created.
    """

    # TaskName.UPDATE_VCF_DIMENSIONS.value is string, with is needed to not have to handle enums in the migrations

    task_history_entry = TaskHistoryDB(
        task_id=None, task_name=TaskName.UPDATE_VCF_DIMENSIONS.value, user_id=user_id, project_id=project_id
    )
    db.add(task_history_entry)
    await db.commit()
    await db.refresh(task_history_entry)
    return task_history_entry.id


async def get_active_dimensions_contenders(
    db: AsyncSession,
    project_id: int,
    provisional_entry_ttl_seconds: int = 120,
    celerytaskmeta_entry_gap_ttl_seconds: int = 3600,  # 1 hour in the queue without being picked up by a worker
) -> list[int]:
    """
    Get active dimensions-update contenders for a project.

    provisional_entry_cutoff is the time between: the task submission arrives to the API and is written to the db as a
    provisional entry (task_id is None), AND the update of the entry with the celery task_id after async_apply has sucessfully enqued the task.

    celery_taskmeta_entry_gap_cutoff is the time between: the celery task_id is written to the db after async_apply has sucessfully enqued the task,
    AND the celery_taskmeta entry is written to the db with the status (i.e. a worker has picked up the task and started).
    """
    now = datetime.now(timezone.utc)
    provisional_entry_cutoff = now - timedelta(seconds=provisional_entry_ttl_seconds)
    celery_taskmeta_entry_gap_cutoff = now - timedelta(seconds=celerytaskmeta_entry_gap_ttl_seconds)

    stmt = (
        select(TaskHistoryDB.id)
        .outerjoin(CeleryTaskMeta, CeleryTaskMeta.task_id == TaskHistoryDB.task_id)
        .where(
            TaskHistoryDB.project_id == project_id,
            TaskHistoryDB.task_name == TaskName.UPDATE_VCF_DIMENSIONS.value,
            or_(
                and_(
                    TaskHistoryDB.task_id.is_(None),
                    TaskHistoryDB.created_at >= provisional_entry_cutoff,
                ),
                and_(
                    TaskHistoryDB.task_id.is_not(None),
                    CeleryTaskMeta.status.in_(ACTIVE_CELERY_STATUSES),
                ),
                and_(
                    TaskHistoryDB.task_id.is_not(None),
                    CeleryTaskMeta.task_id.is_(None),
                    TaskHistoryDB.created_at >= celery_taskmeta_entry_gap_cutoff,
                ),
            ),
        )
        .order_by(TaskHistoryDB.id.asc())
    )

    result = await db.execute(stmt)
    return result.scalars().all()


async def delete_dimensions_provisional_reservation(db: AsyncSession, job_id: int) -> None:
    """
    Deletes a provisional dimensions reservation entry. For instance if it has lost a concurrent race.
    """
    stmt = delete(TaskHistoryDB).where(
        and_(
            TaskHistoryDB.id == job_id,
            TaskHistoryDB.task_name == TaskName.UPDATE_VCF_DIMENSIONS.value,
            TaskHistoryDB.task_id.is_(None),
        )
    )
    await db.execute(stmt)
    await db.commit()
