"""
Unit and integration tests for poll_task_until_final_state_reached retry behaviour.

Fast tests use stamina's testing mode (no real backoff waits) and mock
make_authenticated_request at the module level.

Slow / very-slow tests let stamina run with real backoff delays to verify that
the retry loop actually converges over time.  They require no external
infrastructure — only time.

  pytest                          → fast tests only
  pytest --run-slow               → also runs the ~15 s backoff test
  pytest --run-very-slow          → also runs the ~2 min long-task simulation
"""

from unittest.mock import MagicMock, patch

import pytest
import stamina

from divbase_cli.cli_commands.query_cli import poll_task_until_final_state_reached

DIVBASE_URL = "http://localhost:8001/api"


def _make_task_response(task_id: int, status: str) -> MagicMock:
    """Build a mock HTTP response as returned by make_authenticated_request for a task-history lookup."""
    mock = MagicMock()
    mock.json.return_value = [{"id": task_id, "status": status}]
    return mock


@pytest.fixture
def stamina_no_wait():
    """Disable stamina backoff delays and cap attempts for unit tests that are intended to be fast/not test the actually backoff delays."""
    with stamina.set_testing(testing=True, attempts=15):
        yield


@pytest.mark.usefixtures("stamina_no_wait")
class TestPollingCasesThatAreFast:
    @patch("divbase_cli.cli_commands.query_cli.make_authenticated_request")
    def test_poll_returns_success_when_task_immediately_succeeds(self, mock_make_authenticated_request):
        """Test that poll_task_until_final_state_reached returns 'SUCCESS' on the first poll when the task is already done."""
        mock_make_authenticated_request.return_value = _make_task_response(task_id=7, status="SUCCESS")

        result = poll_task_until_final_state_reached(divbase_url=DIVBASE_URL, task_id=7)

        assert result == "SUCCESS"
        mock_make_authenticated_request.assert_called_once()

    @patch("divbase_cli.cli_commands.query_cli.make_authenticated_request")
    def test_poll_returns_failure_when_task_immediately_fails(self, mock_make_authenticated_request):
        """Test that poll_task_until_final_state_reached returns 'FAILURE' on the first poll when the task already failed."""
        mock_make_authenticated_request.return_value = _make_task_response(task_id=7, status="FAILURE")

        result = poll_task_until_final_state_reached(divbase_url=DIVBASE_URL, task_id=7)

        assert result == "FAILURE"
        mock_make_authenticated_request.assert_called_once()

    @patch("divbase_cli.cli_commands.query_cli.make_authenticated_request")
    def test_poll_retries_through_pending_state_until_success(self, mock_make_authenticated_request):
        """Test that poll_task_until_final_state_reached retries on PENDING responses and returns SUCCESS once the task completes."""
        mock_make_authenticated_request.side_effect = [
            _make_task_response(task_id=7, status="PENDING"),
            _make_task_response(task_id=7, status="PENDING"),
            _make_task_response(task_id=7, status="SUCCESS"),
        ]

        result = poll_task_until_final_state_reached(divbase_url=DIVBASE_URL, task_id=7)

        assert result == "SUCCESS"
        assert mock_make_authenticated_request.call_count == 3

    @patch("divbase_cli.cli_commands.query_cli.make_authenticated_request")
    def test_poll_retries_through_started_state_until_failure(self, mock_make_authenticated_request):
        """Test that poll_task_until_final_state_reached retries on STARTED responses and returns FAILURE when the task ultimately fails."""
        mock_make_authenticated_request.side_effect = [
            _make_task_response(task_id=7, status="STARTED"),
            _make_task_response(task_id=7, status="FAILURE"),
        ]

        result = poll_task_until_final_state_reached(divbase_url=DIVBASE_URL, task_id=7)

        assert result == "FAILURE"
        assert mock_make_authenticated_request.call_count == 2


class TestPollingCasesThatAreSlow:
    @pytest.mark.skipif("not config.getoption('--run-slow')", reason="Only run when --run-slow is given")
    @patch("divbase_cli.cli_commands.query_cli.make_authenticated_request")
    def test_poll_retries_with_real_backoff_until_task_succeeds(self, mock_make_authenticated_request):
        """
        Test that poll_task_until_final_state_reached converges to SUCCESS through real stamina backoff delays.

        Tests four PENDING responses before reaching SUCCESS state. The backoff delays are approximately
        1 + 2 + 4 + 8 ≈ 15 seconds (based on the retries config of poll_task_until_final_state_reached).

        """
        mock_make_authenticated_request.side_effect = [
            _make_task_response(task_id=7, status="PENDING"),
            _make_task_response(task_id=7, status="PENDING"),
            _make_task_response(task_id=7, status="PENDING"),
            _make_task_response(task_id=7, status="PENDING"),
            _make_task_response(task_id=7, status="SUCCESS"),
        ]

        result = poll_task_until_final_state_reached(divbase_url=DIVBASE_URL, task_id=7)

        assert result == "SUCCESS"
        assert mock_make_authenticated_request.call_count == 5


class TestPollingCasesThatAreVerySlow:
    @pytest.mark.skipif("not config.getoption('--run-very-slow')", reason="Only run when --run-very-slow is given")
    @patch("divbase_cli.cli_commands.query_cli.make_authenticated_request")
    def test_poll_handles_long_running_task_with_many_retries(self, mock_make_authenticated_request):
        """
        Test that poll_task_until_final_state_reached handles a long-running task that stays in PENDING
        for many polling cycles before eventually succeeding. Designed to reach the max exponential backoff
        delay of 60 seconds.

        Tests seven PENDING responses before reaching SUCCESS state. The backoff delays are approximately
        1 + 2 + 4 + 8 + 16 + 32 + 60 ≈ 123 seconds (based on the retries config of poll_task_until_final_state_reached).
        """
        pending_responses = [_make_task_response(task_id=7, status="PENDING") for _ in range(7)]
        mock_make_authenticated_request.side_effect = pending_responses + [
            _make_task_response(task_id=7, status="SUCCESS")
        ]

        result = poll_task_until_final_state_reached(divbase_url=DIVBASE_URL, task_id=7)

        assert result == "SUCCESS"
        assert mock_make_authenticated_request.call_count == 8
