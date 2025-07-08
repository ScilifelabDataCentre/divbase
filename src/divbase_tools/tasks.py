import dataclasses
import logging
import os
from pathlib import Path

from celery import Celery

from divbase_tools.cli_commands.query_cli import DEFAULT_METADATA_TSV
from divbase_tools.queries import BCFToolsInput, BcftoolsQueryManager, run_sidecar_metadata_query
from divbase_tools.s3_client import S3FileManager, create_s3_file_manager

BROKER_URL = os.environ.get("CELERY_BROKER_URL", "pyamqp://guest@localhost//")
RESULT_BACKEND = os.environ.get("CELERY_RESULT_BACKEND", "redis://localhost:6379/0")

app = Celery("divbase_tools", broker=BROKER_URL, backend=RESULT_BACKEND)

# Redis-specific config
app.conf.update(
    result_expires=3600,
    task_track_started=True,
    task_serializer="json",
    accept_content=["json"],
    result_serializer="json",
)

logger = logging.getLogger(__name__)


@app.task(name="tasks.sample_metadata_query", tags=["quick"])
def sample_metadata_query_task(tsv_filter: str, bucket_name: str) -> dict:
    """
    Run a sample metadata query task as a Celery task.
    """
    s3_file_manager = create_s3_file_manager(url="http://minio:9000")

    download_sample_metadata(bucket_name=bucket_name, s3_file_manager=s3_file_manager)

    metadata_result = run_sidecar_metadata_query(
        file=DEFAULT_METADATA_TSV,
        filter_string=tsv_filter,
    )
    logger.info(f"Type of metadata_result: {type(metadata_result)}")
    # celery serializes the return value, hence conversion to dict.
    return dataclasses.asdict(metadata_result)


@app.task(name="tasks.bcftools_query", tags=["slow"])
def bcftools_full_query_task(tsv_filter: str, command: str, bucket_name: str, user_name: str = "Default User"):
    """
    Run a full bcftools query command as a Celery task, with sample metadata filtering run first.

    TODO - how MinIO URL is handled needs a lot of thought here.
    Unlike rest of API, might not be in same "local network"
    """
    task_id = bcftools_full_query_task.request.id
    logger.info(f"Starting bcftools_full_query_task with Celery, task ID: {task_id}")

    s3_file_manager = create_s3_file_manager(url="http://minio:9000")

    download_sample_metadata(bucket_name=bucket_name, s3_file_manager=s3_file_manager)

    metadata_result = run_sidecar_metadata_query(
        file=DEFAULT_METADATA_TSV,
        filter_string=tsv_filter,
    )

    _ = download_vcf_files(
        files_to_download=metadata_result.filenames,
        bucket_name=bucket_name,
        s3_file_manager=s3_file_manager,
    )

    bcftools_inputs = dataclasses.asdict(
        BCFToolsInput(
            sample_and_filename_subset=metadata_result.sample_and_filename_subset,
            sampleIDs=metadata_result.sampleIDs,
            filenames=metadata_result.filenames,
        )
    )

    try:
        output_file = BcftoolsQueryManager().execute_pipe(command, bcftools_inputs)
    except Exception as e:
        logger.error(f"Error in bcftools task: {str(e)}")
        return {"status": "error", "error": str(e), "task_id": task_id}

    if bucket_name:
        upload_results_file(output_file=Path(output_file), bucket_name=bucket_name, s3_file_manager=s3_file_manager)

    return {"status": "completed", "output_file": f"{output_file}", "submitter": "celery"}


def download_sample_metadata(bucket_name: str, s3_file_manager: S3FileManager) -> Path:
    """
    Download the metadata file from the specified S3 bucket.
    """
    return s3_file_manager.download_files(
        objects={DEFAULT_METADATA_TSV.name: None},
        download_dir=Path.cwd(),
        bucket_name=bucket_name,
    )[0]


def download_vcf_files(files_to_download: list[str], bucket_name: str, s3_file_manager: S3FileManager) -> list[Path]:
    """
    Fetch input VCF files for bcftools run from the s3 bucket.
    """
    objects = {file_name: None for file_name in files_to_download}
    return s3_file_manager.download_files(
        objects=objects,
        download_dir=Path.cwd(),
        bucket_name=bucket_name,
    )


def upload_results_file(output_file: Path, bucket_name: str, s3_file_manager: S3FileManager) -> None:
    """
    Upon completion of the task, upload the results file to the specified bucket.
    """
    _ = s3_file_manager.upload_files(
        to_upload={output_file.name: output_file},
        bucket_name=bucket_name,
    )
