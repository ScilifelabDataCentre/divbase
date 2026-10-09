"""
Unit tests for the time limits of the worker tasks.
"""

from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest
from celery.exceptions import SoftTimeLimitExceeded

from divbase_api.worker import tasks
from divbase_api.worker.tasks import bcftools_pipe_task, sample_metadata_query_task, update_vcf_dimensions_task
from divbase_lib.exceptions import TaskUserError


def test_dimensions_task_soft_time_limit_gives_max_run_time_error(monkeypatch, tmp_path):
    """The soft time limit should give the user a clear max run time error."""
    s3_file_manager = MagicMock()
    s3_file_manager.list_files.return_value = ["file.vcf.gz"]
    s3_file_manager.latest_version_of_all_files.return_value = {"file.vcf.gz": "v1"}
    monkeypatch.setattr(tasks, "_create_s3_file_manager", lambda: s3_file_manager)
    monkeypatch.setattr(tasks, "SyncSessionLocal", MagicMock())
    monkeypatch.setattr(tasks, "get_vcf_metadata_by_project", lambda **kwargs: SimpleNamespace(vcf_files=[]))
    monkeypatch.setattr(tasks, "get_skipped_vcfs_by_project_worker", lambda **kwargs: {})
    monkeypatch.setattr(tasks, "_remove_stale_dimensions_db_entries", lambda **kwargs: [])
    monkeypatch.setattr(tasks, "_download_vcf_files", lambda **kwargs: {"file.vcf.gz": tmp_path / "file.vcf.gz"})

    def raise_soft_time_limit(self, vcf_path):
        raise SoftTimeLimitExceeded()

    monkeypatch.setattr(tasks.VCFDimensionCalculator, "calculate_dimensions", raise_soft_time_limit)

    with pytest.raises(TaskUserError) as excinfo:
        update_vcf_dimensions_task(bucket_name="bucket", project_id=1, project_name="project", user_id=1)

    error_msg = str(excinfo.value)
    assert "exceeded the maximum run time of" in error_msg
    assert "unexpected error" not in error_msg


def test_bcftools_task_soft_time_limit_gives_max_run_time_error(monkeypatch, tmp_path):
    """The soft time limit should give the user a clear max run time error, and the user's log file is still uploaded."""
    monkeypatch.chdir(tmp_path)  # the task writes its user log file to the current directory
    monkeypatch.setattr(tasks, "_create_s3_file_manager", lambda: MagicMock())
    monkeypatch.setattr(tasks, "SyncSessionLocal", MagicMock())
    monkeypatch.setattr(tasks, "_get_bcftools_version", lambda: "test")

    def raise_soft_time_limit(**kwargs):
        raise SoftTimeLimitExceeded()

    monkeypatch.setattr(tasks, "get_vcf_metadata_by_project", raise_soft_time_limit)

    upload_log_file = MagicMock()
    monkeypatch.setattr(tasks, "_upload_log_file", upload_log_file)

    with pytest.raises(TaskUserError, match="exceeded the maximum run time of"):
        bcftools_pipe_task(
            tsv_filter=None,
            metadata_tsv_name=None,
            command="view",
            bucket_name="bucket",
            project_id=1,
            project_name="project",
            user_id=1,
            job_id=1,
        )

    upload_log_file.assert_called_once()


def test_sample_metadata_query_task_soft_time_limit_gives_max_run_time_error(monkeypatch):
    monkeypatch.setattr(tasks, "_create_s3_file_manager", lambda: MagicMock())

    def raise_soft_time_limit(**kwargs):
        raise SoftTimeLimitExceeded()

    monkeypatch.setattr(tasks, "_download_sample_metadata_for_task", raise_soft_time_limit)

    with pytest.raises(TaskUserError, match="exceeded the maximum run time of"):
        sample_metadata_query_task(
            tsv_filter="Area:North",
            metadata_tsv_name="metadata.tsv",
            bucket_name="bucket",
            project_id=1,
            project_name="project",
            user_id=1,
        )
