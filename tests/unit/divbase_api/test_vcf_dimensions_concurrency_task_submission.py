"""
Unit tests for VCF dimensions API route for the race condition gates to handle task submission concurrency for the same DivBase project.
"""

import asyncio
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from types import SimpleNamespace
from unittest.mock import AsyncMock, MagicMock, patch

from divbase_api.models.projects import ProjectRoles
from divbase_api.routes import vcf_dimensions
from divbase_lib.api_schemas.vcf_dimensions import DimensionsUpdateSubmitResult


@dataclass
class ProjectTest:
    """Dimensions endpoint expects a project object with id, name, and bucket_name attributes, so define a simple dataclass for testing purposes."""

    id: int
    name: str
    bucket_name: str


@dataclass
class UserTest:
    """Dimensions endpoint expects a user object with an id attribute, so define a simple dataclass for testing purposes."""

    id: int


def _run_update_endpoint(db: MagicMock) -> DimensionsUpdateSubmitResult:
    """Helper function that runs the update endpoint with the given mocked db and fixed project/user info. Returns the endpoint result."""
    project = ProjectTest(id=7, name="test-project", bucket_name="test-bucket")
    user = UserTest(id=99)
    with (
        ThreadPoolExecutor(max_workers=1) as pool
    ):  # Use separate thread per test call to avoid event loop collisions/errors when running the full test stack with pytest -s -v
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
    Test that if an active dimensions update job already exists, the endpoint returns that job id and does not create a new provisional entry or enqueues a new task.

    For the first active_dimensions_contenders_since_before gate in the endpoint.
    """
    db = MagicMock()

    with (
        patch.object(
            vcf_dimensions,
            "get_active_dimensions_contenders",
            new=AsyncMock(
                return_value=[42]
            ),  # Simulate hitting the first get_active_dimensions_contenders gate when there already is a job with ID 42
        ) as mock_get_active,
        patch.object(
            vcf_dimensions,
            "create_provisional_dimensions_reservation",
            new=AsyncMock(),
        ) as mock_create_reservation,
        patch.object(vcf_dimensions.update_vcf_dimensions_task, "apply_async") as mock_apply_async,
    ):
        result = _run_update_endpoint(db=db)

    assert result == DimensionsUpdateSubmitResult(job_id=42, outcome="existing")  # Return the existing active job ID 42
    mock_get_active.assert_awaited_once_with(
        db, 7
    )  # There was already an existing job for the project, so get_active_dimensions_contenders should have been called once (first gate)
    mock_create_reservation.assert_not_awaited()  # Provisional reservation should not be created since existing active job was found
    mock_apply_async.assert_not_called()  # No new task should be enqueued since existing active job was found


def test_update_vcf_dimensions_endpoint_loser_deletes_provisional_entry_and_reuses_winner_job():
    """
    Test that if a contender finds an active job but that job is not the contender itself, it deletes its provisional reservation and returns the existing active job as the winner.

    For the second active_dimensions_contenders_since_before gate in the endpoint.
    """
    db = MagicMock()

    with (
        patch.object(
            vcf_dimensions,
            "get_active_dimensions_contenders",
            new=AsyncMock(
                side_effect=[[], [10, 11]]
            ),  # Simulate no active jobs at first get_active_dimensions_contenders gate, then simulate existing active job with ID 10 and contender's provisional job ID 12 at second get_active_dimensions_contenders gate
        ) as mock_get_active,
        patch.object(
            vcf_dimensions,
            "create_provisional_dimensions_reservation",
            new=AsyncMock(return_value=11),
        ) as mock_create_reservation,  # Simulate that the current job creates a provisional reservation with job ID 11
        patch.object(
            vcf_dimensions,
            "delete_dimensions_provisional_reservation",
            new=AsyncMock(),  # Simulate that the current job deletes its provisional reservation
        ) as mock_delete_reservation,
        patch.object(vcf_dimensions.update_vcf_dimensions_task, "apply_async") as mock_apply_async,
    ):
        result = _run_update_endpoint(db=db)

    assert result == DimensionsUpdateSubmitResult(
        job_id=10, outcome="existing"
    )  # Return job ID 10 as the existing active job (the winner), not job ID 11 (the contender)
    assert mock_get_active.await_count == 2
    mock_create_reservation.assert_awaited_once_with(db=db, user_id=99, project_id=7)  #
    mock_delete_reservation.assert_awaited_once_with(
        db=db, job_id=11
    )  # The other job ID 10 is the winner, so job ID 11 should be deleted from the table
    mock_apply_async.assert_not_called()


def test_update_vcf_dimensions_endpoint_winner_enqueues_and_updates_task_history():
    """Test that if a contender does not find any active job other than itself, it enqueues a new task and updates the task history entry with the Celery task ID."""
    db = MagicMock()

    with (
        patch.object(
            vcf_dimensions,
            "get_active_dimensions_contenders",
            new=AsyncMock(side_effect=[[], [33]]),
        ) as mock_get_active,
        patch.object(
            vcf_dimensions,
            "create_provisional_dimensions_reservation",
            new=AsyncMock(return_value=33),
        ) as mock_create_reservation,
        patch.object(
            vcf_dimensions,
            "delete_dimensions_provisional_reservation",
            new=AsyncMock(),
        ) as mock_delete_reservation,
        patch.object(
            vcf_dimensions,
            "update_task_history_entry_with_celery_task_id",
            new=AsyncMock(),
        ) as mock_update_task_history,
        patch.object(
            vcf_dimensions.update_vcf_dimensions_task,
            "apply_async",
            return_value=SimpleNamespace(id="celery-33"),
        ) as mock_apply_async,
    ):
        result = _run_update_endpoint(db=db)

    assert result == DimensionsUpdateSubmitResult(job_id=33, outcome="new")
    assert mock_get_active.await_count == 2
    mock_create_reservation.assert_awaited_once_with(db=db, user_id=99, project_id=7)
    mock_delete_reservation.assert_not_awaited()
    mock_apply_async.assert_called_once_with(
        kwargs={
            "bucket_name": "test-bucket",
            "project_id": 7,
            "project_name": "test-project",
            "user_id": 99,
        }
    )
    mock_update_task_history.assert_awaited_once_with(db=db, job_id=33, task_id="celery-33")
