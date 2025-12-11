"""
Definitions of periodic cron tasks for Celery Beat.
"""

import logging
import os
from datetime import datetime, timedelta, timezone

from celery.schedules import crontab
from sqlalchemy import delete

from divbase_api.models.task_history import CeleryTaskMeta, TaskHistoryDB, TaskStatus
from divbase_api.worker.tasks import app
from divbase_api.worker.worker_db import SyncSessionLocal

logger = logging.getLogger(__name__)


# TODO decide if these should be defined in the env var or implemented in another way. right now they fallback on the default values
TASK_RETENTION_DAYS = int(os.environ.get("TASK_RETENTION_DAYS", "30"))
STUCK_PENDING_STATUS_HOURS = int(os.environ.get("STUCK_PENDING_STATUS_HOURS", "24"))
STUCK_STARTED_STATUS_HOURS = int(os.environ.get("STUCK_STARTED_STATUS_HOURS", "48"))

# TODO NEW QUERY LOGIC NEEDED TO GET PENDING STATUS FROM HERE


@app.task(name="cron_tasks.cleanup_old_task_history")
def cleanup_old_task_history_task(retention_days: int = TASK_RETENTION_DAYS):
    """
    Periodic task to clean up old task history entries from both TaskHistoryDB and CeleryTaskMeta.
    Runs daily to remove entries older than retention_days.
    """
    try:
        cutoff_date = datetime.now(timezone.utc) - timedelta(days=retention_days)

        with SyncSessionLocal() as db:
            deleted_celery_task_meta = db.execute(
                delete(CeleryTaskMeta).where(CeleryTaskMeta.date_done < cutoff_date)
            ).rowcount
            deleted_task_history = db.execute(
                delete(TaskHistoryDB).where(TaskHistoryDB.created_at < cutoff_date)
            ).rowcount

            db.commit()

            logger.info(
                f"Cleaned up {deleted_celery_task_meta} entries from CeleryTaskMeta and "
                f"{deleted_task_history} entries from TaskHistoryDB older than {retention_days} days "
                f"(cutoff: {cutoff_date.isoformat()})"
            )

            return {
                "status": "completed",
                "number_of_celery_meta_deleted": deleted_celery_task_meta,
                "number_of_task_history_deleted": deleted_task_history,
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

    Tasks stuck in PENDING for more than stuck_pending_hours are likely orphaned
    (worker never picked them up or queue issue).

    Tasks stuck in STARTED for more than stuck_started_hours are likely from
    worker crashes where the task started but never completed.
    """
    try:
        now = datetime.now(timezone.utc)
        pending_cutoff = now - timedelta(hours=stuck_pending_hours)
        started_cutoff = now - timedelta(hours=stuck_started_hours)

        with SyncSessionLocal() as db:
            stuck_pending = db.execute(
                delete(TaskHistoryDB)
                .where(TaskHistoryDB.status == TaskStatus.PENDING)
                .where(TaskHistoryDB.created_at < pending_cutoff)
                .returning(TaskHistoryDB.task_id)
            ).fetchall()

            stuck_pending_ids = [row[0] for row in stuck_pending]

            stuck_started = db.execute(
                delete(TaskHistoryDB)
                .where(TaskHistoryDB.status == TaskStatus.STARTED)
                .where(TaskHistoryDB.started_at < started_cutoff)
                .returning(TaskHistoryDB.task_id)
            ).fetchall()

            stuck_started_ids = [row[0] for row in stuck_started]

            if stuck_pending_ids or stuck_started_ids:
                all_stuck_ids = stuck_pending_ids + stuck_started_ids
                db.execute(delete(CeleryTaskMeta).where(CeleryTaskMeta.task_id.in_(all_stuck_ids)))

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
}
