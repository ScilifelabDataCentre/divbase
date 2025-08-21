import dataclasses
import logging
import os
from pathlib import Path

from celery import Celery

from divbase_tools.queries import BCFToolsInput, BcftoolsQueryManager, run_sidecar_metadata_query
from divbase_tools.s3_client import S3FileManager, create_s3_file_manager

logger = logging.getLogger(__name__)

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


def dynamic_router(name, args, kwargs, options, task=None, **kw):
    """
    Function to dynamically route tasks based on their name and arguments.
    The input arguments of the function are taken from the Celery docs.

    Advanced logic can be implemented here to route tasks to different queues
    based on the task names or arguments.

    For now, just use a simple routing based on substrings in the task name.

    To use static routing instead, use the following the call instead of the dynamic_router function:
        app.conf.task_routes = {
        "tasks.simulate_quick_task": {"queue": "quick"},
        "tasks.simulate_long_task": {"queue": "long"},
        "tasks.bcftools_pipe": {"queue": "long"},
        }
    """
    if "quick" in name:
        return {"queue": "quick"}
    if "long" in name:
        return {"queue": "long"}
    if name == "tasks.sample_metadata_query":
        return {"queue": "quick"}
    if name == "tasks.bcftools_query":
        return {"queue": "long"}
    return {"queue": "celery"}


app.conf.task_routes = (dynamic_router,)


@app.task(name="tasks.sample_metadata_query", tags=["quick"])
def sample_metadata_query_task(tsv_filter: str, metadata_tsv_name: str, bucket_name: str) -> dict:
    """
    Run a sample metadata query task as a Celery task.
    """
    s3_file_manager = create_s3_file_manager(url="http://minio:9000")

    metadata_path = download_sample_metadata(
        metadata_tsv_name=metadata_tsv_name, bucket_name=bucket_name, s3_file_manager=s3_file_manager
    )

    metadata_result = run_sidecar_metadata_query(
        file=metadata_path,
        filter_string=tsv_filter,
    )
    # celery serializes the return value, hence conversion to dict.
    return dataclasses.asdict(metadata_result)


@app.task(name="tasks.bcftools_query", tags=["slow"])
def bcftools_pipe_task(
    tsv_filter: str, command: str, metadata_tsv_name: str, bucket_name: str, user_name: str = "Default User"
):
    """
    Run a full bcftools query command as a Celery task, with sample metadata filtering run first.

    TODO - how MinIO URL is handled needs a lot of thought here.
    Unlike rest of API, might not be in same "local network"
    """
    task_id = bcftools_pipe_task.request.id
    logger.info(f"Starting bcftools_pipe_task with Celery, task ID: {task_id}")

    s3_file_manager = create_s3_file_manager(url="http://minio:9000")

    metadata_path = download_sample_metadata(
        metadata_tsv_name=metadata_tsv_name, bucket_name=bucket_name, s3_file_manager=s3_file_manager
    )

    metadata_result = run_sidecar_metadata_query(
        file=metadata_path,
        filter_string=tsv_filter,
    )

    _ = download_vcf_files(
        files_to_download=metadata_result.unique_filenames,
        bucket_name=bucket_name,
        s3_file_manager=s3_file_manager,
    )

    bcftools_inputs = dataclasses.asdict(
        BCFToolsInput(
            sample_and_filename_subset=metadata_result.sample_and_filename_subset,
            sampleIDs=metadata_result.unique_sample_ids,
            filenames=metadata_result.unique_filenames,
        )
    )

    try:
        output_file = BcftoolsQueryManager().execute_pipe(command, bcftools_inputs, task_id)
    except Exception as e:
        logger.error(f"Error in bcftools task: {str(e)}")
        return {"status": "error", "error": str(e), "task_id": task_id}

    upload_results_file(output_file=Path(output_file), bucket_name=bucket_name, s3_file_manager=s3_file_manager)
    delete_job_files_from_worker(
        vcf_paths=metadata_result.unique_filenames, metadata_path=metadata_path, output_file=output_file
    )
    return {"status": "completed", "output_file": output_file, "submitter": user_name}


def download_sample_metadata(metadata_tsv_name: str, bucket_name: str, s3_file_manager: S3FileManager) -> Path:
    """
    Download the metadata file from the specified S3 bucket.
    """
    return s3_file_manager.download_files(
        objects={metadata_tsv_name: None},
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


def delete_job_files_from_worker(vcf_paths: list[Path], metadata_path: Path, output_file: Path) -> None:
    """
    After uploading results to bucket, delete job files from the worker.
    """
    for vcf_path in vcf_paths:
        try:
            os.remove(vcf_path)
            logger.info(f"deleted {vcf_path}")
        except Exception as e:
            logger.warning(f"Could not delete input VCF file from worker {vcf_path}: {e}")
    try:
        os.remove(metadata_path)
        logger.info(f"deleted {metadata_path}")
    except Exception as e:
        logger.warning(f"Could not delete metadata file from worker {metadata_path}: {e}")
    try:
        os.remove(output_file)
        logger.info(f"deleted {output_file}")
    except Exception as e:
        logger.warning(f"Could not delete output file from worker {output_file}: {e}")
