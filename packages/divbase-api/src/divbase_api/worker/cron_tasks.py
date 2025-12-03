"""
Definitions of periodic cron tasks for Celery Beat.
"""

import logging
import os
from datetime import datetime, timedelta, timezone

from celery.schedules import crontab
from sqlalchemy import delete

from divbase_api.models.task_history import CeleryTaskMeta, TaskHistoryDB
from divbase_api.worker.tasks import app
from divbase_api.worker.worker_db import SyncSessionLocal

logger = logging.getLogger(__name__)


TASK_RETENTION_DAYS = int(os.environ.get("TASK_RETENTION_DAYS", "30"))


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
                "celery_meta_deleted": deleted_celery_task_meta,
                "task_history_deleted": deleted_task_history,
                "cutoff_date": cutoff_date.isoformat(),
                "retention_days": retention_days,
            }

    except Exception as e:
        logger.error(f"Failed to cleanup old task history: {e}")
        raise


app.conf.beat_schedule = {
    "cleanup-old-tasks-daily": {
        "task": "cron_tasks.cleanup_old_task_history",
        "schedule": crontab(hour=5, minute=0),  # Run daily at 5 AM. Don't set to 2 AM or 3 AM due to daylight saving
        "kwargs": {"retention_days": TASK_RETENTION_DAYS},
    },
}
