"""
Definitions of periodic cron tasks for Celery Beat.
"""

import logging
import os
from datetime import datetime, timedelta, timezone

from celery.schedules import crontab
from sqlalchemy import delete, select, text

from divbase_api.models.project_versions import ProjectVersionDB
from divbase_api.models.projects import ProjectDB
from divbase_api.models.revoked_tokens import RevokedTokenDB
from divbase_api.models.task_history import CeleryTaskMeta, TaskHistoryDB, TaskStartedAtDB
from divbase_api.services.s3_client import create_s3_file_manager
from divbase_api.worker.tasks import S3_ENDPOINT_URL, app
from divbase_api.worker.worker_db import SyncSessionLocal
from divbase_lib.divbase_constants import QUERY_RESULTS_FILE_PREFIX

logger = logging.getLogger(__name__)


# TODO decide if these should be defined in the env var or implemented in another way. right now they fallback on the default values
TASK_RETENTION_DAYS = int(os.environ.get("TASK_RETENTION_DAYS", "30"))
STUCK_PENDING_STATUS_HOURS = int(os.environ.get("STUCK_PENDING_STATUS_HOURS", "168"))  # 168 h = 7 days
STUCK_STARTED_STATUS_HOURS = int(os.environ.get("STUCK_STARTED_STATUS_HOURS", "168"))  # 168 h = 7 days

SOFT_DELETED_FILES_RETENTION_DAYS = 30
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


@app.task(name="cron_tasks.update_storage_usage_metrics")
def update_storage_usage_metrics():
    """
    Periodic task to update storage usage metrics for all projects.

    This task runs daily and calculates the total storage used by each project, including all versions and files.
    """
    s3_file_manager = create_s3_file_manager(url=S3_ENDPOINT_URL)
    with SyncSessionLocal() as db:
        stmt = select(ProjectDB).where(ProjectDB.is_active == True)  # noqa: E712
        projects = db.execute(stmt).scalars().all()

        for project in projects:
            storage_used_bytes = s3_file_manager.get_bucket_usage_bytes(bucket_name=project.bucket_name)
            project.storage_used_bytes = storage_used_bytes
            db.add(project)

        db.commit()

    return {
        "status": "completed",
        "number_of_projects_updated": len(projects),
    }


@app.task(name="cron_tasks.hard_delete_expired_soft_deleted_objects")
def hard_delete_expired_soft_deleted_objects():
    """
    Periodic task to hard delete any soft-deleted objects from each project's S3 bucket after given retention period.

    All versions of a file that has been soft deleted for more than 30 days are hard deleted:
    This rule has 1 exception, if a version of a file is recorded in a user defined project version, then it is not hard deleted.

    In those cases, a delete marker is re-added to the file. Otherwise, the undeleted files would show up
    as the latest version of the file and be queryable/usable, which would not be expected behaviour.

    NOTE: Results files are not handled here, they are handled in a seperate task that also deletes the associated job entries in the db tables.
    """
    cutoff_date = datetime.now(timezone.utc) - timedelta(days=SOFT_DELETED_FILES_RETENTION_DAYS)
    s3_file_manager = create_s3_file_manager(url=S3_ENDPOINT_URL)

    per_project_delete_count = {}
    # Mapping of protected file versions per project ID, to avoid deleting
    protected_project_files: dict[int, dict[str, set[str]]] = {}
    with SyncSessionLocal() as db:
        projects = db.execute(statement=select(ProjectDB)).scalars().all()
        for project in projects:
            per_project_delete_count[project.name] = 0

            protected_versions_map: dict[str, set[str]] = {}
            stmt = select(ProjectVersionDB).where(ProjectVersionDB.project_id == project.id)
            project_versions = db.execute(stmt).scalars().all()

            for version in project_versions:
                if not version.files:
                    continue
                for filename, version_id in version.files.items():
                    if filename not in protected_versions_map:
                        protected_versions_map[filename] = set()
                    protected_versions_map[filename].add(version_id)

            protected_project_files[project.id] = protected_versions_map

            # look for objects where the current file version is a delete marker older than the cutoff
            candidate_objects_to_purge = s3_file_manager.get_expired_soft_deleted_objects(
                bucket_name=project.bucket_name,
                cutoff_date=cutoff_date,
                prefix_exclude=QUERY_RESULTS_FILE_PREFIX,  # handled by diff cron task.
            )

            for object_name in candidate_objects_to_purge:
                protected_ids = protected_versions_map.get(object_name, set())

                # We need to list all versions of this specific key to determine what to keep/delete
                # easier to use the s3_client directly
                versions_resp = s3_file_manager.s3_client.list_object_versions(
                    Bucket=project.bucket_name, Prefix=object_name
                )

                # Delete both file versions and delete markers of the object.
                objects_to_delete = []
                for version in versions_resp.get("Versions", []):
                    if version["Key"] == object_name and version["VersionId"] not in protected_ids:
                        objects_to_delete.append({"Key": object_name, "VersionId": version["VersionId"]})

                for marker in versions_resp.get("DeleteMarkers", []):
                    if marker["Key"] == object_name:
                        objects_to_delete.append({"Key": object_name, "VersionId": marker["VersionId"]})

                if objects_to_delete:
                    s3_file_manager.hard_delete_specific_object_versions(
                        objects=objects_to_delete,
                        bucket_name=project.bucket_name,
                    )

                    per_project_delete_count[project.name] += len(objects_to_delete)

            # We take advantage of a special behaviour of S3, that you can delete objects that don't exist.
            # This ensures that if a protected version now becomes the latest, it wont be used in queries/downloadable etc..
            s3_file_manager.soft_delete_objects(
                objects=list(candidate_objects_to_purge),
                bucket_name=project.bucket_name,
            )

    return {
        "status": "completed",
        "objects_hard_deleted_per_project": per_project_delete_count,
        "cutoff_date": cutoff_date.isoformat(),
        "retention_days": SOFT_DELETED_FILES_RETENTION_DAYS,
    }


# NOTE! If you add a new task here, make sure it starts with "cron_tasks"
# Don't set to 2 AM or 3 AM due to daylight saving.
# Timezone for job schedule is CET (defined in app in tasks.py).
app.conf.beat_schedule = {
    "cleanup-old-tasks-daily": {
        "task": "cron_tasks.cleanup_old_task_history",
        "schedule": crontab(hour=5, minute=0),  # Run daily at 5 AM CET
        "kwargs": {"retention_days": TASK_RETENTION_DAYS},
    },
    "cleanup-stuck-tasks-daily": {
        "task": "cron_tasks.cleanup_stuck_tasks",
        "schedule": crontab(hour=5, minute=15),  # Run daily at 5:15 AM CET
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
    "update-storage-usage-metrics-daily": {
        "task": "cron_tasks.update_storage_usage_metrics",
        "schedule": crontab(hour=5, minute=30),  # Run daily at 5:30 AM CET
    },
    "hard-delete-expired-soft-deleted-objects-weekly": {
        "task": "cron_tasks.hard_delete_expired_soft_deleted_objects",
        "schedule": crontab(day_of_week=1, hour=5, minute=35),  # Run every Monday at 5:35 AM CET
    },
}
