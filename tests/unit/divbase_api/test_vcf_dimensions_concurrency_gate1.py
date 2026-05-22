"""
Unit tests for the Gate 1 concurrency check in the VCF dimensions update endpoint.
Gate 1: if an active dimensions update job already exists for the project, return it instead of enqueueing a new one.
"""

import asyncio
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from types import SimpleNamespace
from unittest.mock import AsyncMock, MagicMock, patch

from divbase_api.models.projects import ProjectRoles
from divbase_api.routes import vcf_dimensions
from divbase_api.worker.task_names import TaskName
from divbase_lib.api_schemas.vcf_dimensions import DimensionsUpdateSubmitResult


@dataclass
class ProjectTest:
    """Minimal project object matching what the endpoint unpacks from project_and_user_and_role."""

    id: int
    name: str
    bucket_name: str


@dataclass
class UserTest:
    """Minimal user object matching what the endpoint unpacks from project_and_user_and_role."""

    id: int
    is_admin: bool = False


def _run_update_endpoint(db: MagicMock) -> DimensionsUpdateSubmitResult:
    """Run the update endpoint in a separate thread to avoid event loop collisions when running the full test suite."""
    project = ProjectTest(id=7, name="test-project", bucket_name="test-bucket")
    user = UserTest(id=99)
    with ThreadPoolExecutor(max_workers=1) as pool:
        return pool.submit(
            lambda: asyncio.run(
                vcf_dimensions.update_vcf_dimensions_endpoint(
                    project_name=project.name,
                    project_and_user_and_role=(project, user, ProjectRoles.EDIT),
                    db=db,
                )
            )
        ).result(timeout=5.0)


def test_update_vcf_dimensions_endpoint_reuses_existing_active_job():
    """
    Gate 1: if an active dimensions update job already exists for the project, the endpoint returns
    that job id immediately without enqueueing a new task or writing to task history.
    """
    db = MagicMock()

    with (
        patch.object(vcf_dimensions, "check_queue_closed_for_new_tasks", new=AsyncMock()),
        patch.object(
            vcf_dimensions,
            "get_active_dimensions_contenders",
            new=AsyncMock(return_value=[42]),
        ) as mock_get_active,
        patch.object(vcf_dimensions, "create_task_history_entry", new=AsyncMock()) as mock_create_history,
        patch.object(vcf_dimensions.update_vcf_dimensions_task, "apply_async") as mock_apply_async,
    ):
        result = _run_update_endpoint(db=db)

    assert result == DimensionsUpdateSubmitResult(job_id=42, outcome="existing")
    mock_get_active.assert_awaited_once_with(db=db, project_id=7)
    mock_create_history.assert_not_awaited()
    mock_apply_async.assert_not_called()


def test_update_vcf_dimensions_endpoint_enqueues_new_job_when_no_active_job_exists():
    """
    Gate 1: if no active dimensions update job exists, the endpoint enqueues a new task,
    records it in task history with task_name set, and returns outcome="new".
    """
    db = MagicMock()

    with (
        patch.object(vcf_dimensions, "check_queue_closed_for_new_tasks", new=AsyncMock()),
        patch.object(
            vcf_dimensions,
            "get_active_dimensions_contenders",
            new=AsyncMock(return_value=[]),
        ) as mock_get_active,
        patch.object(
            vcf_dimensions.update_vcf_dimensions_task,
            "apply_async",
            return_value=SimpleNamespace(id="celery-uuid-33"),
        ) as mock_apply_async,
        patch.object(
            vcf_dimensions,
            "create_task_history_entry",
            new=AsyncMock(return_value=33),
        ) as mock_create_history,
    ):
        result = _run_update_endpoint(db=db)

    assert result == DimensionsUpdateSubmitResult(job_id=33, outcome="new")
    mock_get_active.assert_awaited_once_with(db=db, project_id=7)
    mock_apply_async.assert_called_once_with(
        kwargs={
            "bucket_name": "test-bucket",
            "project_id": 7,
            "project_name": "test-project",
            "user_id": 99,
        }
    )
    mock_create_history.assert_awaited_once_with(
        user_id=99,
        project_id=7,
        task_id="celery-uuid-33",
        task_name=TaskName.UPDATE_VCF_DIMENSIONS.value,
        db=db,
    )
