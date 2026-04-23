"""
Integration tests for task-history behavior that assert backend internals.
"""

import time

from sqlalchemy import select
from typer.testing import CliRunner

from divbase_api.models.task_history import CeleryTaskMeta, TaskHistoryDB
from divbase_api.worker.worker_db import SyncSessionLocal
from divbase_cli.divbase_cli import app

runner = CliRunner()


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
                    .where(TaskHistoryDB.id == task_id)
                )
                result = db.execute(stmt).scalar_one_or_none()

                if result is not None:
                    break

                time.sleep(retry_delay)

            assert result is not None, f"Task {task_id} not found in results backend after {max_retries} retries"
            assert result in ["PENDING", "STARTED", "SUCCESS", "FAILURE"]
