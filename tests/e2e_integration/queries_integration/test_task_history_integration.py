"""
Integration tests for task-history behavior that assert backend internals.
"""

import time
import uuid
from csv import reader
from io import StringIO

from sqlalchemy import delete, select
from typer.testing import CliRunner

from divbase_api.models.task_history import CeleryTaskMeta, TaskHistoryDB
from divbase_api.models.users import UserDB
from divbase_api.worker.tasks import BCFTOOLS_QUERY_TASK_NAME, SAMPLE_METADATA_QUERY_TASK_NAME
from divbase_api.worker.worker_db import SyncSessionLocal
from divbase_cli.divbase_cli import app

runner = CliRunner()


def _parse_tsv_rows(stdout: str) -> list[list[str]]:
    """Parse task-history --tsv output into rows, skipping the header"""
    return list(reader(StringIO(stdout), delimiter="\t"))[1:]


def test_get_task_status_by_task_id_uses_results_backend(
    CONSTANTS, logged_in_edit_user_with_existing_config, run_update_dimensions, project_map
):
    """
    Verify that user task IDs returned from CLI submissions can be resolved to Celery states
    via the PostgreSQL results backend.
    """
    project_name = CONSTANTS["QUERY_PROJECT"]
    project_id = project_map[project_name]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    user_id = 1
    run_update_dimensions(bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id)

    tsv_filter = "Area:West of Ireland,Northern Portugal;"
    arg_command = "view -r 21:15000000-25000000"
    command = f"query vcf --tsv-filter '{tsv_filter}' --command '{arg_command}' --project {project_name} "

    first_task_result = runner.invoke(app, command)
    assert first_task_result.exit_code == 0
    first_task_id = first_task_result.stdout.strip().split()[-1]

    second_task_result = runner.invoke(app, command)
    assert second_task_result.exit_code == 0
    second_task_id = second_task_result.stdout.strip().split()[-1]

    max_retries = 10
    retry_delay = 0.5

    with SyncSessionLocal() as db:
        for task_id in [first_task_id, second_task_id]:
            result = None
            for _ in range(max_retries):
                stmt = (
                    select(CeleryTaskMeta.status)
                    .join(TaskHistoryDB, CeleryTaskMeta.task_id == TaskHistoryDB.task_id)
                    .where(TaskHistoryDB.id == int(task_id))
                )
                result = db.execute(stmt).scalar_one_or_none()

                if result is not None:
                    break

                time.sleep(retry_delay)

            assert result is not None, f"Task {task_id} not found in results backend after {max_retries} retries"
            assert result in ["PENDING", "STARTED", "SUCCESS", "FAILURE"]


def test_queued_tasks_without_celery_meta_are_shown_as_queuing(
    CONSTANTS, logged_in_edit_user_with_existing_config, project_map
):
    """
    A task waiting in the queue only has a TaskHistoryDB entry. The CeleryTaskMeta and TaskStartedAtDB entries
    are only created once a worker picks up the task.
    Theses tasks should be shown as QUEUING.

    To simulate this state we add TaskHistoryDB entries directly (otherwise worker would pick them up).
    """
    project_name = CONSTANTS["QUERY_PROJECT"]
    project_id = project_map[project_name]
    user_email = CONSTANTS["TEST_USERS"]["edit user"]["email"]

    with SyncSessionLocal() as db:
        user_id = db.execute(select(UserDB.id).where(UserDB.email == user_email)).scalar_one()
        queued_tasks = [
            # initial state when created by query tsv task
            TaskHistoryDB(
                task_id=f"{uuid.uuid4()}",
                user_id=user_id,
                project_id=project_id,
                task_name=SAMPLE_METADATA_QUERY_TASK_NAME,
            ),
            TaskHistoryDB(
                task_id=f"{uuid.uuid4()}",
                user_id=user_id,
                project_id=project_id,
                task_name=SAMPLE_METADATA_QUERY_TASK_NAME,
            ),
            # initial state when created by query vcf task
            TaskHistoryDB(task_id=None, user_id=user_id, project_id=project_id, task_name=BCFTOOLS_QUERY_TASK_NAME),
        ]
        db.add_all(queued_tasks)
        db.commit()
        queued_job_ids = {task.id for task in queued_tasks}  # divbase job id (aka rolling ints)

    try:
        # validate divbase-cli task-history id
        for job_id in queued_job_ids:
            result = runner.invoke(app, f"task-history id {job_id} --tsv")
            assert result.exit_code == 0, result.output
            rows = _parse_tsv_rows(result.stdout)
            assert len(rows) == 1
            assert rows[0][1] == str(job_id)
            assert rows[0][2] == "QUEUING"

        # validate divbase-cli task-history user --project
        result = runner.invoke(app, f"task-history user --project {project_name} --tsv")
        assert result.exit_code == 0, result.output
        rows = _parse_tsv_rows(result.stdout)
        states_by_job_id = {int(row[1]): row[2] for row in rows}
        assert queued_job_ids <= states_by_job_id.keys()
        for job_id in queued_job_ids:
            assert states_by_job_id[job_id] == "QUEUING"

        # validate divbase-cli task-history user
        result = runner.invoke(app, "task-history user --tsv")
        assert result.exit_code == 0, result.output
        rows = _parse_tsv_rows(result.stdout)
        states_by_job_id = {int(row[1]): row[2] for row in rows}
        assert queued_job_ids <= states_by_job_id.keys()
        for job_id in queued_job_ids:
            assert states_by_job_id[job_id] == "QUEUING"
    finally:
        with SyncSessionLocal() as db:
            db.execute(delete(TaskHistoryDB).where(TaskHistoryDB.id.in_(queued_job_ids)))
            db.commit()
