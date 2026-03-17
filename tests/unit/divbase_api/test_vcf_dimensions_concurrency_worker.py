"""
Unit tests for worker-side dimensions race condition gates to handle task execution concurrency for the same DivBase project.
"""

from unittest.mock import MagicMock, patch

from divbase_api.worker.crud_dimensions import resolve_dimensions_winner_for_worker
from divbase_api.worker.tasks import _check_if_concurrent_dimensions_update_task_exist
from divbase_lib.api_schemas.vcf_dimensions import DimensionUpdateTaskResult


def test_resolve_dimensions_winner_for_worker_defaults_to_current_when_no_winner_found():
    """
    Test that if the worker does not find any active contenders other than the current job,
    it treats the current job as the winner by default to fail open and avoid blocking dimensions updates.
    """
    db = MagicMock()
    mock_result_1 = MagicMock()
    mock_result_1.scalar_one_or_none.return_value = (
        9  # Simulate finding only one active contender with the same job ID as the current job
    )
    mock_result_2 = MagicMock()
    mock_result_2.scalar_one_or_none.return_value = (
        None  # Simulate not finding any active contender with a different job ID than the current job
    )
    db.execute.side_effect = [mock_result_1, mock_result_2]

    is_winner, current_job_id, winner_job_id = resolve_dimensions_winner_for_worker(
        db=db,
        project_id=2,
        task_id="celery-task-2",
    )

    assert is_winner is True
    assert current_job_id == 9
    assert winner_job_id == 9
    assert db.execute.call_count == 2


def test_resolve_dimensions_winner_for_worker_returns_non_winner_when_other_job_wins():
    """
    Test that if the worker finds an active contender with a different job ID than the current job,
    it treats the current job as the non-winner and returns the winning job ID.
    """
    db = MagicMock()
    mock_result_1 = MagicMock()
    mock_result_1.scalar_one_or_none.return_value = 9  # Simulate finding an active contender with a different job ID (9) than the current job (which will be simulated as 3 in the second query)
    mock_result_2 = MagicMock()
    mock_result_2.scalar_one_or_none.return_value = 3  # Simulate finding an active contender with a different job ID (3) than the current job (which is simulated as 9 in the first query)
    db.execute.side_effect = [mock_result_1, mock_result_2]

    is_winner, current_job_id, winner_job_id = resolve_dimensions_winner_for_worker(
        db=db,
        project_id=2,
        task_id="celery-task-3",
    )

    assert is_winner is False
    assert current_job_id == 9
    assert winner_job_id == 3
    assert db.execute.call_count == 2


def test_check_if_concurrent_dimensions_update_task_exist_returns_none_for_winner():
    """
    Test that if the worker determines that the current task is the winner among any active contenders for dimensions update,
    the function returns None to indicate that the task can proceed with updating dimensions.
    """
    with (
        patch("divbase_api.worker.tasks.SyncSessionLocal") as mock_session_local,
        patch(
            "divbase_api.worker.tasks.resolve_dimensions_winner_for_worker",
            return_value=(
                True,
                11,
                11,
            ),  # Simulate that the current task is the winner with job ID 11 among active contenders
        ) as mock_resolve_winner,
    ):
        mock_session_local.return_value.__enter__.return_value = MagicMock()
        result = _check_if_concurrent_dimensions_update_task_exist(project_id=7, task_id="celery-task-4")

    assert result is None
    mock_resolve_winner.assert_called_once()


def test_check_if_concurrent_dimensions_update_task_exist_returns_skipped_duplicate_result():
    """
    Test that if the worker determines that the current task is not the winner among active contenders for dimensions update and finds a winning job ID,
    it returns a DimensionUpdateTaskResult indicating that the task is a skipped duplicate along with the winning job ID.
    """
    with (
        patch("divbase_api.worker.tasks.SyncSessionLocal") as mock_session_local,
        patch(
            "divbase_api.worker.tasks.resolve_dimensions_winner_for_worker",
            return_value=(
                False,
                12,
                5,
            ),  # Simulate that the current task is not the winner, with its own job ID being 12 and the winning job ID being 5
        ) as mock_resolve_winner,
    ):
        mock_session_local.return_value.__enter__.return_value = MagicMock()
        result = _check_if_concurrent_dimensions_update_task_exist(project_id=7, task_id="celery-task-5")

    assert result is not None
    assert isinstance(result, DimensionUpdateTaskResult)
    assert result.status == "skipped_duplicate"
    assert result.duplicate_of_job_id == 5
    assert result.message == "Skipped duplicate dimensions update; active job id: 5"
    mock_resolve_winner.assert_called_once()


def test_resolve_current_task_dimensions_winner_for_worker_when_task_history_missing():
    """
    Test that that the edge case where the task history entry for the current task is missing when the worker tries to resolve a winner for concurrent dimensions update tasks,
    it defaults to treating the current task as the winner to avoid blocking dimensions updates ("fail open" strategy).
    """
    db = MagicMock()
    mock_result = MagicMock()
    mock_result.scalar_one_or_none.return_value = None  # Simulate missing entry in TaskHistoryDB for the current task
    db.execute.side_effect = [mock_result]

    is_winner, current_job_id, winner_job_id = resolve_dimensions_winner_for_worker(
        db=db,
        project_id=1,
        task_id="celery-task-1",
    )

    assert is_winner is True
    assert current_job_id is None
    assert winner_job_id is None
    assert db.execute.call_count == 1
