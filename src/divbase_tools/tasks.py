import dataclasses
import logging
import os
import re
from pathlib import Path

import yaml
from celery import Celery

from divbase_tools.exceptions import NoVCFFilesFoundError
from divbase_tools.queries import BCFToolsInput, BcftoolsQueryManager, run_sidecar_metadata_query
from divbase_tools.s3_client import S3FileManager, create_s3_file_manager
from divbase_tools.vcf_dimension_indexing import VCFDimensionIndexManager

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

    if "view -r" in command:
        files_to_download = check_for_unnecessary_files_for_region_query(
            bucket_name=bucket_name,
            s3_file_manager=s3_file_manager,
            command=command,
            files_to_download=metadata_result.unique_filenames,
        )

        sample_and_filename_subset = [
            entry for entry in metadata_result.sample_and_filename_subset if entry["Filename"] in files_to_download
        ]
    else:
        sample_and_filename_subset = metadata_result.sample_and_filename_subset

    _ = download_vcf_files(
        files_to_download=files_to_download,
        bucket_name=bucket_name,
        s3_file_manager=s3_file_manager,
    )

    bcftools_inputs = dataclasses.asdict(
        BCFToolsInput(
            sample_and_filename_subset=sample_and_filename_subset,
            sampleIDs=metadata_result.unique_sample_ids,
            filenames=files_to_download,
        )
    )

    try:
        output_file = BcftoolsQueryManager().execute_pipe(command, bcftools_inputs)
    except Exception as e:
        logger.error(f"Error in bcftools task: {str(e)}")
        return {"status": "error", "error": str(e), "task_id": task_id}

    upload_results_file(output_file=Path(output_file), bucket_name=bucket_name, s3_file_manager=s3_file_manager)
    return {"status": "completed", "output_file": output_file, "submitter": user_name}


@app.task(name="tasks.update_vcf_dimensions_task")
def update_vcf_dimensions_task(bucket_name: str, user_name: str = "Default User"):
    task_id = update_vcf_dimensions_task.request.id
    s3_file_manager = create_s3_file_manager(url="http://minio:9000")

    all_files = s3_file_manager.list_files(bucket_name=bucket_name)
    vcf_files = [file for file in all_files if file.endswith(".vcf") or file.endswith(".vcf.gz")]

    if not vcf_files:
        raise NoVCFFilesFoundError(
            f"VCF dimensions file could not be generated since no VCF files were found in the project: {bucket_name}. Please upload at least one VCF file and run this command again."
        )

    manager = VCFDimensionIndexManager(bucket_name=bucket_name, s3_file_manager=s3_file_manager)
    already_indexed_vcfs = manager.get_indexed_filenames()

    non_indexed_vcfs = [file for file in vcf_files if file not in already_indexed_vcfs]

    _ = download_vcf_files(
        files_to_download=non_indexed_vcfs,
        bucket_name=bucket_name,
        s3_file_manager=s3_file_manager,
    )

    files_indexed_by_this_job = []
    for file in non_indexed_vcfs:
        try:
            manager.update_dimension_entry(vcf_filename=file)
            files_indexed_by_this_job.append(file)
        except Exception as e:
            logger.error(f"Error in dimensions indexing task: {str(e)}")
            return {"status": "error", "error": str(e), "task_id": task_id}
    return {
        "status": "completed",
        "submitter": user_name,
        "VCF files that were added to dimensions file by this job": files_indexed_by_this_job,
    }

    # TODO delete downloaded files from worker upon fininshing


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


def check_for_unnecessary_files_for_region_query(
    bucket_name: str, s3_file_manager: S3FileManager, command: str, files_to_download: list[str]
) -> dict:
    """
    If the 'view -r' query is present, read the .vcf_dimensions.yaml file and check if the specified scaffolds are available in the VCF files.
    If a file in files_to_download (as identified by sidecar sample metadata query) does not contain any of the specified scaffolds, it will be skipped.
    This will save time and resources by avoiding unnecessary transfer of irrelevant files between the bucket and the worker.

    NOTE! There is a risk that users mistype...

    TODO ensure that the errors raised here are raised before submitting the task to celery.
    """

    try:
        dimensions_index = read_vcf_dimensions_file(bucket_name=bucket_name, s3_file_manager=s3_file_manager)
    except FileNotFoundError:
        return files_to_download
        # TODO for now, this makes it so that if the dimensions file is not present, the job continues by skipping the check of unnecessary files. There should at least be a message that tells this to the user

    scaffolds = []
    files_to_download_updated = []
    matches = re.findall(r"view\s+-r\s+([^\s;]+)", command)
    for match in matches:
        scaffolds.extend([region.strip() for region in match.split(",") if region.strip()])

    for file in files_to_download:
        record = None
        for rec in dimensions_index.get("dimensions", []):
            if rec.get("filename") == file:
                record = rec
                break
        if not record:
            print(f"File '{file}' is not indexed in the VCF dimensions file.")
            continue

        record_scaffolds = set(record.get("dimensions", {}).get("scaffolds", []))
        for scaffold in scaffolds:
            if scaffold in record_scaffolds:
                print(f"Scaffold '{scaffold}' is present in file '{file}'.")
                if file not in files_to_download_updated:
                    files_to_download_updated.append(file)

    if files_to_download_updated == []:
        raise ValueError(
            "Based on the 'view -r' query and the VCF scaffolds indexed in DivBase, there are no VCF files in the project that fulfills the query. Please try another -r query with scaffolds/chromosomes that are present in the VCF files."
        )
        # TODO this stops jobs from executing. Which might not be intuitive. Such subsests would be able to run in bftools but the result would be an merged VCF with no variants (header only). Should this be a warning instead that will be returned in the worker results?

    return files_to_download_updated

    # TODO when there is a command implemented that reads the .vcf_dimensions.yaml file and list all scaffolds in the project, add that command as a help to this error message.


def read_vcf_dimensions_file(bucket_name: str, s3_file_manager: S3FileManager) -> dict:
    try:
        content = s3_file_manager.download_s3_file_to_str(key=".vcf_dimensions.yaml", bucket_name=bucket_name)
    except Exception as exc:
        raise FileNotFoundError(f".vcf_dimensions.yaml not found in bucket '{bucket_name}'.") from exc
    # TODO add a hint to the error telling how dimension indexing can be done.
    # TODO use the get_dimensions_info instance method of the VCFDimensionIndexManager class instead of reading file diretly? don't want to download the file though
    return yaml.safe_load(content)
