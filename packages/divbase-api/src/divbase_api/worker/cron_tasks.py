"""
Definitions of periodic cron tasks for Celery Beat.
"""

import logging
import os
from datetime import datetime, timedelta, timezone

from celery.schedules import crontab
from sqlalchemy import delete, text

from divbase_api.models.project_versions import ProjectVersionDB
from divbase_api.models.revoked_tokens import RevokedTokenDB
from divbase_api.models.task_history import CeleryTaskMeta, TaskHistoryDB, TaskStartedAtDB
from divbase_api.worker.tasks import app
from divbase_api.worker.worker_db import SyncSessionLocal

logger = logging.getLogger(__name__)


# TODO decide if these should be defined in the env var or implemented in another way. right now they fallback on the default values
TASK_RETENTION_DAYS = int(os.environ.get("TASK_RETENTION_DAYS", "30"))
STUCK_PENDING_STATUS_HOURS = int(os.environ.get("STUCK_PENDING_STATUS_HOURS", "168"))  # 168 h = 7 days
STUCK_STARTED_STATUS_HOURS = int(os.environ.get("STUCK_STARTED_STATUS_HOURS", "168"))  # 168 h = 7 days

SOFT_DELETED_PROJECT_VERSION_RETENTION_DAYS = 30
REVOKED_TOKEN_MAX_AGE_DAYS = 7


@app.task(name="cron_tasks.cleanup_old_task_history")
def cleanup_old_task_history_task(retention_days: int = TASK_RETENTION_DAYS):
    """
    Periodic task to clean up old task history entries from both TaskHistoryDB and CeleryTaskMeta.
    Runs daily to remove entries older than retention_days.
    """
    try:
        cutoff_date = datetime.now(timezone.utc) - timedelta(days=retention_days)

        with SyncSessionLocal() as db:
            old_task_ids = [
                row[0]
                for row in db.execute(
                    text("SELECT task_id FROM task_history WHERE created_at < :cutoff_date"),
                    {"cutoff_date": cutoff_date},
                ).fetchall()
            ]

            deleted_celery_task_meta = db.execute(
                delete(CeleryTaskMeta).where(CeleryTaskMeta.task_id.in_(old_task_ids))
            ).rowcount
            deleted_task_history = db.execute(
                delete(TaskHistoryDB).where(TaskHistoryDB.task_id.in_(old_task_ids))
            ).rowcount
            deleted_started_at = db.execute(
                delete(TaskStartedAtDB).where(TaskStartedAtDB.task_id.in_(old_task_ids))
            ).rowcount

            db.commit()

            logger.info(
                f"Cleaned up {deleted_celery_task_meta} entries from CeleryTaskMeta, "
                f"{deleted_task_history} from TaskHistoryDB, and "
                f"{deleted_started_at} from TaskStartedAtDB older than {retention_days} days "
                f"(cutoff: {cutoff_date.isoformat()})"
            )

            return {
                "status": "completed",
                "number_of_celery_meta_deleted": deleted_celery_task_meta,
                "number_of_task_history_deleted": deleted_task_history,
                "number_of_started_at_deleted": deleted_started_at,
                "cutoff_date": cutoff_date.isoformat(),
                "retention_days": retention_days,
            }

    except Exception as e:
        logger.error(f"Failed to cleanup old task history: {e}")
        raise


@app.task(name="cron_tasks.cleanup_stuck_tasks")
def cleanup_stuck_tasks_task(
    stuck_pending_hours: int = STUCK_PENDING_STATUS_HOURS,
    stuck_started_hours: int = STUCK_STARTED_STATUS_HOURS,
):
    """
    Periodic task to clean up tasks stuck in PENDING or STARTED status from
    TaskHistoryDB and CeleryTaskMeta.

    Joins TaskHistoryDB and CeleryTaskMeta, finds tasks stuck in PENDING or STARTED
    for longer than the configured threshold, and deletes them from both tables.

    Database lookup is based on celery task ID, since that is a UUID that is written to
    all tables.
    """
    try:
        now = datetime.now(timezone.utc)
        pending_cutoff = now - timedelta(hours=stuck_pending_hours)
        started_cutoff = now - timedelta(hours=stuck_started_hours)

        with SyncSessionLocal() as db:
            # Find stuck PENDING tasks
            stuck_pending_ids = [
                row[0]
                for row in db.execute(
                    text("""
                    SELECT th.task_id
                    FROM task_history th
                    JOIN celery_taskmeta cm ON th.task_id = cm.task_id
                    WHERE cm.status = 'PENDING'
                    AND th.created_at < :pending_cutoff
                    """),
                    {"pending_cutoff": pending_cutoff},
                ).fetchall()
            ]

            # Find stuck STARTED tasks
            stuck_started_ids = [
                row[0]
                for row in db.execute(
                    text("""
                    SELECT th.task_id
                    FROM task_history th
                    JOIN celery_taskmeta cm ON th.task_id = cm.task_id
                    JOIN task_started_at tsa ON th.task_id = tsa.task_id
                    WHERE cm.status = 'STARTED'
                    AND tsa.started_at < :started_cutoff
                    """),
                    {"started_cutoff": started_cutoff},
                ).fetchall()
            ]

            all_stuck_ids = stuck_pending_ids + stuck_started_ids

            if all_stuck_ids:
                db.execute(delete(TaskHistoryDB).where(TaskHistoryDB.task_id.in_(all_stuck_ids)))
                db.execute(delete(CeleryTaskMeta).where(CeleryTaskMeta.task_id.in_(all_stuck_ids)))
                db.execute(delete(TaskStartedAtDB).where(TaskStartedAtDB.task_id.in_(all_stuck_ids)))
                db.commit()

            logger.info(
                f"Cleaned up {len(stuck_pending_ids)} PENDING tasks older than {stuck_pending_hours}h "
                f"and {len(stuck_started_ids)} STARTED tasks older than {stuck_started_hours}h"
            )

            return {
                "status": "completed",
                "number_of_stuck_pending_deleted": len(stuck_pending_ids),
                "number_of_stuck_started_deleted": len(stuck_started_ids),
                "pending_cutoff": pending_cutoff.isoformat(),
                "started_cutoff": started_cutoff.isoformat(),
            }

    except Exception as e:
        logger.error(f"Failed to cleanup stuck tasks: {e}")
        raise


@app.task(name="cron_tasks.cleanup_old_revoked_tokens")
def cleanup_old_revoked_tokens():
    """
    Periodic task to clean up old revoked token entries.
    (These tokens will all have expired by this timepoint anyway.)
    """
    cutoff_date = datetime.now(timezone.utc) - timedelta(days=REVOKED_TOKEN_MAX_AGE_DAYS)
    stmt = delete(RevokedTokenDB).where(RevokedTokenDB.revoked_at < cutoff_date)
    with SyncSessionLocal() as db:
        deleted_count = db.execute(stmt).rowcount
        db.commit()
    return {
        "status": "completed",
        "number_of_revoked_tokens_deleted": deleted_count,
        "cutoff_date": cutoff_date.isoformat(),
        "max_revoked_token_age_days": REVOKED_TOKEN_MAX_AGE_DAYS,
    }


@app.task(name="cron_tasks.cleanup_soft_deleted_project_versions")
def cleanup_soft_deleted_project_versions():
    """
    Periodic task to hard delete any soft-deleted project versions older than the retention period.
    """
    cutoff_date = datetime.now(timezone.utc) - timedelta(days=SOFT_DELETED_PROJECT_VERSION_RETENTION_DAYS)
    stmt = delete(ProjectVersionDB)
    stmt = stmt.where(ProjectVersionDB.is_deleted == True)  # noqa: E712
    stmt = stmt.where(ProjectVersionDB.date_deleted < cutoff_date)

    with SyncSessionLocal() as db:
        deleted_count = db.execute(stmt).rowcount
        db.commit()
    return {
        "status": "completed",
        "number_of_project_versions_hard_deleted": deleted_count,
        "cutoff_date": cutoff_date.isoformat(),
        "soft_delete_retention_period_days": SOFT_DELETED_PROJECT_VERSION_RETENTION_DAYS,
    }


# NOTE! If you add a new task here, make sure it starts with "cron_tasks"
app.conf.beat_schedule = {
    "cleanup-old-tasks-daily": {
        "task": "cron_tasks.cleanup_old_task_history",
        "schedule": crontab(
            hour=5, minute=0
        ),  # Run daily at 5 AM CET (timezone defined in app in tasks.py). Don't set to 2 AM or 3 AM due to daylight saving
        "kwargs": {"retention_days": TASK_RETENTION_DAYS},
    },
    "cleanup-stuck-tasks-daily": {
        "task": "cron_tasks.cleanup_stuck_tasks",
        "schedule": crontab(
            hour=5, minute=15
        ),  # Run daily at 5:15 AM CET (timezone defined in app in tasks.py). Don't set to 2 AM or 3 AM due to daylight saving
        "kwargs": {
            "stuck_pending_hours": STUCK_PENDING_STATUS_HOURS,
            "stuck_started_hours": STUCK_STARTED_STATUS_HOURS,
        },
    },
    "cleanup-old-revoked-daily": {
        "task": "cron_tasks.cleanup_old_revoked_tokens",
        "schedule": crontab(hour=5, minute=20),  # Run daily at 5:20 AM CET
    },
    "cleanup-soft-deleted-project-versions-daily": {
        "task": "cron_tasks.cleanup_soft_deleted_project_versions",
        "schedule": crontab(hour=5, minute=25),  # Run daily at 5:25 AM CET
    },
}
