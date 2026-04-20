"""
Unit tests for task history deserialization service.
"""

import json
from datetime import datetime, timezone

from divbase_api.services.task_history import _deserialize_celery_task_metadata


def test_deserialize_sample_metadata_task_with_legacy_result_keys():
    """Legacy task results with Sample_ID/Filename keys should still deserialize."""
    task = {
        "user_task_id": 42,
        "submitter_email": "user@example.com",
        "status": "SUCCESS",
        "result": json.dumps(
            {
                "sample_and_filename_subset": [
                    {"Sample_ID": "5a_HOM-I7", "Filename": "HOM_20ind_17SNPs_first_10_samples.vcf.gz"},
                    {"Sample_ID": "8_HOM-E57", "Filename": "HOM_20ind_17SNPs_last_10_samples.vcf.gz"},
                ],
                "unique_sample_ids": ["5a_HOM-I7", "8_HOM-E57"],
                "unique_filenames": [
                    "HOM_20ind_17SNPs_first_10_samples.vcf.gz",
                    "HOM_20ind_17SNPs_last_10_samples.vcf.gz",
                ],
                "query_message": "mock",
                "warnings": [],
                "status": "success",
            }
        ),
        "name": "tasks.sample_metadata_query",
        "kwargs": json.dumps(
            {
                "tsv_filter": "Area:West of Ireland",
                "metadata_tsv_name": "sample_metadata.tsv",
                "bucket_name": "mock-bucket",
                "project_id": 1,
                "project_name": "mock-project",
                "user_id": 1,
            }
        ),
        "args": "[]",
        "worker": "celery@worker-1",
        "created_at": datetime(2026, 4, 20, 10, 0, 0, tzinfo=timezone.utc),
        "started_at": datetime(2026, 4, 20, 10, 0, 1, tzinfo=timezone.utc),
        "date_done": datetime(2026, 4, 20, 10, 0, 3, tzinfo=timezone.utc),
    }

    deserialized = _deserialize_celery_task_metadata(task)

    assert deserialized.result is not None
    assert deserialized.result.sample_and_filename_subset[0].sample_id == "5a_HOM-I7"
    assert deserialized.result.sample_and_filename_subset[0].filename == "HOM_20ind_17SNPs_first_10_samples.vcf.gz"
    assert deserialized.result.sample_and_filename_subset[1].sample_id == "8_HOM-E57"
    assert deserialized.result.sample_and_filename_subset[1].filename == "HOM_20ind_17SNPs_last_10_samples.vcf.gz"
