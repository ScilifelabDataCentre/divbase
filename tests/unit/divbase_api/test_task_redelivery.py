"""
Unit tests for failing fast when a long task is redelivered after it was already started (e.g. its worker died).
"""

import pytest

from divbase_api.worker.tasks import bcftools_pipe_task, update_vcf_dimensions_task
from divbase_lib.exceptions import TaskUserError

# task will fail earlier than when any of these values matter, so they don't have to make sense
BCFTOOLS_QUERY_KWARGS = {
    "tsv_filter": None,
    "metadata_tsv_name": None,
    "command": "view",
    "bucket_name": "bucket",
    "project_id": 1,
    "project_name": "project",
    "user_id": 1,
    "job_id": 1,
}
DIMENSIONS_UPDATE_KWARGS = {"bucket_name": "bucket", "project_id": 1, "project_name": "project", "user_id": 1}


@pytest.mark.parametrize(
    "task, task_kwargs",
    [(bcftools_pipe_task, BCFTOOLS_QUERY_KWARGS), (update_vcf_dimensions_task, DIMENSIONS_UPDATE_KWARGS)],
)
def test_long_task_fails_fast_if_task_already_ran(task, task_kwargs):
    """Test that the 2 long-running task types we have fail fast if they have already been run."""
    task.push_request(id="test-task-id", task_already_ran=True)
    try:
        with pytest.raises(TaskUserError, match="Please resubmit the job"):
            task.run(**task_kwargs)
    finally:
        task.pop_request()
