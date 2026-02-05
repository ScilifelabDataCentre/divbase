import dataclasses
import logging
import os
import re
import time
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import psutil
from celery import Celery
from celery.signals import (
    after_task_publish,
    task_prerun,
    worker_process_init,
)

from divbase_api.exceptions import VCFDimensionsEntryMissingError
from divbase_api.models.task_history import TaskHistoryDB, TaskStartedAtDB
from divbase_api.services.queries import BCFToolsInput, BcftoolsQueryManager, run_sidecar_metadata_query
from divbase_api.services.s3_client import S3FileManager, create_s3_file_manager
from divbase_api.worker.crud_dimensions import (
    create_or_update_skipped_vcf,
    create_or_update_vcf_metadata,
    delete_skipped_vcf,
    delete_vcf_metadata,
    get_skipped_vcfs_by_project_worker,
    get_vcf_metadata_by_project,
)
from divbase_api.worker.metrics import (
    ENABLE_WORKER_METRICS_PER_TASK,
    MemoryMonitor,
    TaskMetrics,
    start_metrics_server,
    store_task_metric_in_cache,
    update_prometheus_gauges_from_cache,
)
from divbase_api.worker.vcf_dimension_indexing import (
    VCFDimensionCalculator,
)
from divbase_api.worker.worker_db import SyncSessionLocal
from divbase_lib.api_schemas.vcf_dimensions import DimensionUpdateTaskResult
from divbase_lib.exceptions import NoVCFFilesFoundError

logger = logging.getLogger(__name__)


BROKER_URL = os.environ.get("CELERY_BROKER_URL", "pyamqp://guest@localhost//")
RESULT_BACKEND = os.environ.get(
    "CELERY_RESULT_BACKEND", "db+postgresql://divbase_user:badpassword@localhost:5432/divbase_db"
)
S3_ENDPOINT_URL = os.environ.get("S3_ENDPOINT_URL", "http://host.docker.internal:9000")
app = Celery("divbase_worker", broker=BROKER_URL, backend=RESULT_BACKEND)


CELERY_TASKMETA_TABLE_NAME = "celery_taskmeta"
CELERY_GROUPMETA_TABLE_NAME = "celery_groupmeta"

# Celery results backend config
app.conf.update(
    task_track_started=True,
    task_serializer="json",
    accept_content=["json"],
    result_serializer="json",
    result_extended=True,
    timezone="Europe/Stockholm",  # for internal scheduling, e.g. celery beat
    # let celery auto-create db tables (alembic is configured to not manage changes to these tables)
    database_table_names={
        "task": CELERY_TASKMETA_TABLE_NAME,
        "group": CELERY_GROUPMETA_TABLE_NAME,
    },
    result_expires=None,  # disables celery.backend_cleanup since Divbase uses custom cleanup tasks (see cron_tasks.py).
)


@worker_process_init.connect
def init_worker_metrics(**kwargs):
    """
    Start the Prometheus metrics server for the Celery worker. This signal is triggered once per forked worker process.
    With celery prefork concurrency=1, only one worker process will start the metrics server.
    """
    if ENABLE_WORKER_METRICS_PER_TASK:
        start_metrics_server(port=8101)
    else:
        logger.info("Worker metrics collection disabled (ENABLE_WORKER_METRICS_PER_TASK=false)")


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
    if name == "tasks.update_vcf_dimensions_task":
        return {"queue": "long"}  # can take minutes for large VCF files
    return {"queue": "celery"}


app.conf.task_routes = (dynamic_router,)


@after_task_publish.connect
def task_pending_handler(sender=None, headers=None, body=None, **kwargs):
    """
    For cron_tasks: create TaskHistoryDB entry when task is published to broker.

    NOTE! User-submitted tasks are created by the API. This signal handler is for system-submitted tasks.

    """
    task_name = headers["task"]
    if not task_name.startswith("cron_tasks"):
        return

    task_id = headers["id"]

    task_history_entry = TaskHistoryDB(task_id=task_id, user_id=None, project_id=None)
    with SyncSessionLocal() as db:
        db.add(task_history_entry)
        db.commit()


@task_prerun.connect
def handle_task_started(sender=None, task_id=None, **kwargs):
    """
    Signal handler that inserts the started_at timestamp in TaskStartedAtDB
    for the given task_id.
    """
    started_at = datetime.now(timezone.utc)
    if task_id is None:
        raise ValueError("task_id is None in task_prerun signal handler")

    with SyncSessionLocal() as db:
        entry = TaskStartedAtDB(task_id=task_id, started_at=started_at)
        db.add(entry)
        db.commit()


@app.task(name="tasks.sample_metadata_query", tags=["quick"])
def sample_metadata_query_task(
    tsv_filter: str,
    metadata_tsv_name: str,
    bucket_name: str,
    project_id: int,
    project_name: str,
    user_id: int,
) -> dict:
    """Run a sample metadata query task as a Celery task."""
    task_id = sample_metadata_query_task.request.id

    s3_file_manager = create_s3_file_manager(url=S3_ENDPOINT_URL)

    metadata_path = _download_sample_metadata(
        metadata_tsv_name=metadata_tsv_name, bucket_name=bucket_name, s3_file_manager=s3_file_manager
    )

    with SyncSessionLocal() as db:
        vcf_dimensions_data = get_vcf_metadata_by_project(project_id=project_id, db=db)

    if not vcf_dimensions_data.get("vcf_files"):
        raise VCFDimensionsEntryMissingError(project_name=project_name)

    metadata_result = run_sidecar_metadata_query(
        file=metadata_path,
        filter_string=tsv_filter,
        project_id=project_id,
        vcf_dimensions_data=vcf_dimensions_data,
    )

    try:
        os.remove(metadata_path)
        logger.info(f"Deleted metadata file {metadata_path} from worker.")
    except Exception as e:
        logger.warning(f"Could not delete metadata file {metadata_path}: {e}")

    # Convert to dict since celery serializes to JSON when sending back to API layer
    result = metadata_result.model_dump()
    result["status"] = "completed"
    result["task_id"] = task_id

    logger.info(
        f"Metadata query completed: {len(metadata_result.unique_sample_ids)} samples "
        f"mapped to {len(metadata_result.unique_filenames)} VCF files"
    )

    return result


@app.task(name="tasks.bcftools_query", tags=["slow"])
def bcftools_pipe_task(
    tsv_filter: str,
    command: str,
    metadata_tsv_name: str,
    bucket_name: str,
    project_id: int,
    project_name: str,
    user_id: int,
    job_id: int,
):
    """
    Run a full bcftools query command as a Celery task, with sample metadata filtering run first.
    """
    if ENABLE_WORKER_METRICS_PER_TASK:
        task_walltime_start = time.time()
        process = psutil.Process()
        cpu_start = process.cpu_times()
        memory_monitor = MemoryMonitor(process, sample_interval=0.5)
        memory_monitor.start()

    task_id = bcftools_pipe_task.request.id
    logger.info(f"Starting bcftools_pipe_task with Celery, task ID: {task_id}, job ID: {job_id}")

    s3_file_manager = create_s3_file_manager(url=S3_ENDPOINT_URL)

    with SyncSessionLocal() as db:
        vcf_dimensions_data = get_vcf_metadata_by_project(project_id=project_id, db=db)

    if not vcf_dimensions_data.get("vcf_files"):
        raise VCFDimensionsEntryMissingError(project_name=project_name)

    latest_versions_of_bucket_files = s3_file_manager.latest_version_of_all_files(bucket_name=bucket_name)

    metadata_path = _download_sample_metadata(
        metadata_tsv_name=metadata_tsv_name, bucket_name=bucket_name, s3_file_manager=s3_file_manager
    )

    metadata_result = run_sidecar_metadata_query(
        file=metadata_path,
        filter_string=tsv_filter,
        project_id=project_id,
        vcf_dimensions_data=vcf_dimensions_data,
    )

    _check_that_file_versions_match_dimensions_index(
        vcf_dimensions_data=vcf_dimensions_data,
        latest_versions_of_bucket_files=latest_versions_of_bucket_files,
        metadata_result=metadata_result,
    )

    files_to_download = metadata_result.unique_filenames
    sample_and_filename_subset = metadata_result.sample_and_filename_subset

    if "view -r" in command:
        files_to_download = _check_for_unnecessary_files_for_region_query(
            command=command,
            files_to_download=metadata_result.unique_filenames,
            vcf_dimensions_data=vcf_dimensions_data,
        )

        sample_and_filename_subset = [
            entry for entry in metadata_result.sample_and_filename_subset if entry["Filename"] in files_to_download
        ]

    _check_if_samples_can_be_combined_with_bcftools(
        files_to_download=files_to_download,
        vcf_dimensions_data=vcf_dimensions_data,
    )

    # Measure download operation resources
    # Memory: Capture baseline before download, then measure incremental memory increase (delta)
    if ENABLE_WORKER_METRICS_PER_TASK:
        download_cpu_start = process.cpu_times()
        download_mem_baseline = process.memory_info().rss
        download_memory_monitor = MemoryMonitor(process, sample_interval=0.5, baseline_memory=download_mem_baseline)
        download_memory_monitor.start()
        download_start = time.time()

    logger.info("Started downloading VCF files from S3 to worker")

    _ = _download_vcf_files(
        files_to_download=files_to_download,
        bucket_name=bucket_name,
        s3_file_manager=s3_file_manager,
    )

    logger.info("Finished downloading VCF files from S3 to worker")

    if ENABLE_WORKER_METRICS_PER_TASK:
        download_walltime = time.time() - download_start
        download_memory_stats = download_memory_monitor.stop()
        download_cpu_end = process.cpu_times()
        download_cpu_used = (download_cpu_end.user + download_cpu_end.system) - (
            download_cpu_start.user + download_cpu_start.system
        )
        download_mem_peak_delta = download_memory_stats["peak_bytes"]
        download_mem_avg_delta = download_memory_stats["avg_bytes"]
        logger.debug(
            f"VCF download from S3 to worker took (walltime): {download_walltime:.2f}s, CPU: {download_cpu_used:.2f}s"
        )

    bcftools_inputs = dataclasses.asdict(
        BCFToolsInput(
            sample_and_filename_subset=sample_and_filename_subset,
            sampleIDs=metadata_result.unique_sample_ids,
            filenames=files_to_download,
        )
    )

    logger.info("Started bcftools subprocesses")

    # Execute bcftools and get true subprocess metrics
    # Let validation exceptions (BcftoolsPipeEmptyCommandError, BcftoolsPipeUnsupportedCommandError,
    # SidecarInvalidFilterError) propagate to mark task as FAILURE. Otherwise the tasks will incorrectly be marked as SUCCESS.
    output_file, bcftools_metrics = BcftoolsQueryManager().execute_pipe(command, bcftools_inputs, job_id)

    logger.info("Finished bcftools subprocesses")

    bcftools_cpu_used = bcftools_metrics.get("cpu_seconds", 0.0)
    bcftools_mem_peak = bcftools_metrics.get("peak_memory_bytes", 0)
    bcftools_mem_avg = bcftools_metrics.get("avg_memory_bytes", 0)
    bcftools_walltime = bcftools_metrics.get("walltime_seconds", 0.0)

    _upload_results_file(output_file=Path(output_file), bucket_name=bucket_name, s3_file_manager=s3_file_manager)
    _delete_job_files_from_worker(vcf_paths=files_to_download, metadata_path=metadata_path, output_file=output_file)

    if ENABLE_WORKER_METRICS_PER_TASK:
        memory_stats = memory_monitor.stop()
        cpu_end = process.cpu_times()
        python_overhead_cpu = (cpu_end.user + cpu_end.system) - (cpu_start.user + cpu_start.system)
        total_task_cpu = python_overhead_cpu + bcftools_cpu_used
        task_walltime = time.time() - task_walltime_start

        task_metrics = TaskMetrics(
            job_id=job_id,
            task_name=bcftools_pipe_task.name,
            task_cpu_seconds=total_task_cpu,
            task_python_overhead_cpu_seconds=python_overhead_cpu,
            task_memory_peak_bytes=memory_stats["peak_bytes"],
            task_memory_avg_bytes=memory_stats["avg_bytes"],
            task_bcftools_cpu_seconds=bcftools_cpu_used,
            task_bcftools_memory_peak_bytes=bcftools_mem_peak,
            task_bcftools_memory_avg_bytes=bcftools_mem_avg,
            task_bcftools_walltime_seconds=bcftools_walltime,
            task_vcf_download_walltime_seconds=download_walltime,
            task_vcf_download_cpu_seconds=download_cpu_used,
            task_vcf_download_memory_peak_bytes=download_mem_peak_delta,
            task_vcf_download_memory_avg_bytes=download_mem_avg_delta,
            task_walltime_seconds=task_walltime,
        )
        _record_task_metrics(task_metrics)

    return {"status": "completed", "output_file": output_file}


@app.task(name="tasks.update_vcf_dimensions_task")
def update_vcf_dimensions_task(
    bucket_name: str,
    project_id: int,
    project_name: str,
    user_id: int,
) -> dict:
    """
    Update VCF dimensions in the database for the specified bucket.
    """
    task_id = update_vcf_dimensions_task.request.id

    s3_file_manager = create_s3_file_manager(url=S3_ENDPOINT_URL)

    all_files = s3_file_manager.list_files(bucket_name=bucket_name)
    vcf_files = [file for file in all_files if file.endswith(".vcf") or file.endswith(".vcf.gz")]

    if not vcf_files:
        raise NoVCFFilesFoundError(
            f"VCF dimensions file could not be generated since no VCF files were found in the project: {project_name}."
            "Please upload at least one VCF file and run this command again."
        )

    with SyncSessionLocal() as db:
        vcf_dimensions_data = get_vcf_metadata_by_project(project_id=project_id, db=db)

    already_indexed_vcfs = {
        entry["vcf_file_s3_key"]: entry["s3_version_id"] for entry in vcf_dimensions_data.get("vcf_files", [])
    }

    with SyncSessionLocal() as db:
        already_skipped_vcfs = get_skipped_vcfs_by_project_worker(db=db, project_id=project_id)

    latest_versions_of_bucket_files = s3_file_manager.latest_version_of_all_files(bucket_name=bucket_name)

    non_indexed_vcfs = [
        file
        for file in vcf_files
        if not (
            (file in already_skipped_vcfs and already_skipped_vcfs[file] == latest_versions_of_bucket_files.get(file))
            or (
                file in already_indexed_vcfs and already_indexed_vcfs[file] == latest_versions_of_bucket_files.get(file)
            )
        )
    ]

    _ = _download_vcf_files(
        files_to_download=non_indexed_vcfs,
        bucket_name=bucket_name,
        s3_file_manager=s3_file_manager,
    )

    calculator = VCFDimensionCalculator()

    files_indexed_by_this_job = []
    divbase_results_files_skipped_by_this_job = []

    for file in non_indexed_vcfs:
        try:
            vcf_dims = calculator.calculate_dimensions(Path(file))

            if vcf_dims is None:
                skipped_vcf_data = {
                    "vcf_file_s3_key": file,
                    "project_id": project_id,
                    "s3_version_id": latest_versions_of_bucket_files.get(file),
                    "skip_reason": "divbase_generated",
                }

                with SyncSessionLocal() as db:
                    create_or_update_skipped_vcf(db=db, skipped_vcf_data=skipped_vcf_data)

                divbase_results_files_skipped_by_this_job.append(file)
                logger.info(f"Skipping DivBase-generated result file: {file}")
                continue

            vcf_metadata_data = {
                "vcf_file_s3_key": file,
                "project_id": project_id,
                "s3_version_id": latest_versions_of_bucket_files.get(file),
                "samples": vcf_dims["sample_names"],
                "scaffolds": vcf_dims["scaffolds"],
                "variant_count": vcf_dims["variants"],
                "sample_count": vcf_dims["sample_count"],
                "file_size_bytes": Path(file).stat().st_size if Path(file).exists() else 0,
            }

            with SyncSessionLocal() as db:
                create_or_update_vcf_metadata(db=db, vcf_metadata_data=vcf_metadata_data)
            files_indexed_by_this_job.append(file)
            logger.info(f"Indexed VCF metadata for: {file}")

        except Exception as e:
            logger.error(f"Error indexing {file}: {str(e)}")
            return {"status": "error", "error": str(e), "task_id": task_id}

    vcfs_deleted_from_bucket_since_last_indexing = list(set(already_indexed_vcfs) - set(vcf_files))
    if not vcfs_deleted_from_bucket_since_last_indexing:
        vcfs_deleted_from_bucket_since_last_indexing = None
    skipped_deleted_from_bucket = list(set(already_skipped_vcfs) - set(vcf_files))

    # TODO this block could be done in one go for calling db once with a list
    if vcfs_deleted_from_bucket_since_last_indexing:
        for file in vcfs_deleted_from_bucket_since_last_indexing:
            try:
                with SyncSessionLocal() as db:
                    delete_vcf_metadata(db=db, vcf_file_s3_key=file, project_id=project_id)
                logger.info(f"Deleted VCF metadata for removed file: {file}")
            except Exception as e:
                logger.error(f"Failed to delete VCF metadata for {file}: {e}")

    if skipped_deleted_from_bucket:
        for file in skipped_deleted_from_bucket:
            with SyncSessionLocal() as db:
                delete_skipped_vcf(db=db, vcf_file_s3_key=file, project_id=project_id)

    _delete_job_files_from_worker(vcf_paths=non_indexed_vcfs)

    if not files_indexed_by_this_job:
        files_indexed_by_this_job = None
    if not divbase_results_files_skipped_by_this_job:
        divbase_results_files_skipped_by_this_job = None

    result = DimensionUpdateTaskResult(
        status="completed",
        VCF_files_added=files_indexed_by_this_job,
        VCF_files_skipped=divbase_results_files_skipped_by_this_job,
        VCF_files_deleted=vcfs_deleted_from_bucket_since_last_indexing,
    )

    # Convert to dict since celery serializes to JSON when writing to results backend. Pydantic model serialization is not supported by celery
    return result.model_dump()


def _download_sample_metadata(metadata_tsv_name: str, bucket_name: str, s3_file_manager: S3FileManager) -> Path:
    """
    Download the metadata file from the specified S3 bucket.
    """
    return s3_file_manager.download_files(
        objects={metadata_tsv_name: None},
        download_dir=Path.cwd(),
        bucket_name=bucket_name,
    )[0]


def _download_vcf_files(files_to_download: list[str], bucket_name: str, s3_file_manager: S3FileManager) -> list[Path]:
    """
    Fetch input VCF files for bcftools run from the s3 bucket.
    """
    logger.debug(f"Starting download of {len(files_to_download)} VCF file(s) from bucket '{bucket_name}'")

    objects = {file_name: None for file_name in files_to_download}
    downloaded_files = s3_file_manager.download_files(
        objects=objects,
        download_dir=Path.cwd(),
        bucket_name=bucket_name,
    )

    logger.info(f"Downloaded VCF files: {[f.name for f in downloaded_files]}")

    return downloaded_files


def _upload_results_file(output_file: Path, bucket_name: str, s3_file_manager: S3FileManager) -> None:
    """
    Upon completion of the task, upload the results file to the specified bucket.
    """
    _ = s3_file_manager.upload_files(
        to_upload={output_file.name: output_file},
        bucket_name=bucket_name,
    )


def _check_for_unnecessary_files_for_region_query(
    command: str,
    files_to_download: list[str],
    vcf_dimensions_data: dict,
) -> list[str]:
    """
    If the 'view -r' query is present, read the .vcf_dimensions.yaml file and check if the specified scaffolds are available in the VCF files.
    If a file in files_to_download (as identified by sidecar sample metadata query) does not contain any of the specified scaffolds, it will be skipped.
    This will save time and resources by avoiding unnecessary transfer of irrelevant files between the bucket and the worker.

    NOTE! There is a risk that users mistype...

    """
    if not vcf_dimensions_data.get("vcf_files"):
        logger.warning(
            "VCF dimensions data is missing or empty. "
            "All current VCF files will be transferred to the worker without filtering."
        )
        return files_to_download

    scaffolds = []
    files_to_download_updated = []

    matches = re.findall(r"view\s+-r\s+([^\s;]+)", command)
    if not matches:
        return files_to_download

    for match in matches:
        regions = [region.strip() for region in match.split(",") if region.strip()]
        for region in regions:
            scaffold_name = region.split(":")[0]
            scaffolds.append(scaffold_name)

    vcf_lookup = {entry["vcf_file_s3_key"]: entry for entry in vcf_dimensions_data["vcf_files"]}

    for file in files_to_download:
        if file not in vcf_lookup:
            logger.warning(f"File '{file}' is not indexed in the VCF dimensions.")
            continue

        record_scaffolds = set(vcf_lookup[file].get("scaffolds", []))

        for scaffold_name in scaffolds:
            if scaffold_name in record_scaffolds:
                logger.info(f"'view -r' query requires scaffold '{scaffold_name}'. It is present in file '{file}'.")
                if file not in files_to_download_updated:
                    files_to_download_updated.append(file)

    if not files_to_download_updated:
        raise ValueError(
            "Based on the 'view -r' query and the VCF scaffolds indexed in DivBase, "
            "there are no VCF files in the project that fulfill the query.\n"
            "Please try another -r query with scaffolds/chromosomes that are present in the VCF files.\n"
            "To see a list of all unique scaffolds: "
            "'divbase-cli dimensions show --unique-scaffolds --project <PROJECT_NAME>'"
        )

    return files_to_download_updated


def _delete_job_files_from_worker(
    vcf_paths: list[Path] = None, metadata_path: Path = None, output_file: Path = None
) -> None:
    """
    After uploading results to bucket, delete job files from the worker.
    """
    for vcf_path in vcf_paths:
        try:
            os.remove(vcf_path)
            logger.info(f"Deleted {vcf_path} from worker.")
        except Exception as e:
            logger.warning(f"Could not delete input VCF file from worker {vcf_path}: {e}")
    if metadata_path is not None:
        try:
            os.remove(metadata_path)
            logger.info(f"Deleted {metadata_path} from worker.")
        except Exception as e:
            logger.warning(f"Could not delete metadata file from worker {metadata_path}: {e}")
    if output_file is not None:
        try:
            os.remove(output_file)
            logger.info(f"Deleted {output_file} from worker.")
        except Exception as e:
            logger.warning(f"Could not delete output file from worker {output_file}: {e}")


def _check_if_samples_can_be_combined_with_bcftools(
    files_to_download: list[str],
    vcf_dimensions_data: dict,
) -> None:
    """
    Check if samples in VCF files can be combined with bcftools merge/concat.
    Raises ValueError if samples have incompatible overlaps.
    """

    vcf_lookup = {entry["vcf_file_s3_key"]: entry for entry in vcf_dimensions_data["vcf_files"]}

    file_to_samples = {}
    for file in files_to_download:
        if file not in vcf_lookup:
            raise ValueError(f"Sample names not found for file '{file}' in dimensions index.")

        sample_names = vcf_lookup[file].get("samples")
        if not sample_names:
            raise ValueError(f"Sample names not found for file '{file}' in dimensions index.")

        file_to_samples[file] = sample_names

    manager = BcftoolsQueryManager()
    sample_sets = manager._group_vcfs_by_sample_set(file_to_samples)
    logger.debug(f"Sample sets found in the VCF files: {sample_sets}")

    sample_set_overlap_results = _calculate_pairwise_overlap_types_for_sample_sets(sample_sets)
    logger.debug(f"Sample sets overlap type: {sample_set_overlap_results}")

    if (
        sample_set_overlap_results["identical elements, different order"]
        or sample_set_overlap_results["partly overlapping"]
    ):
        msg_lines = []
        if sample_set_overlap_results["identical elements, different order"]:
            msg_lines.append(
                "Sample sets with identical elements but different order:\n"
                + "\n".join(
                    [
                        f"{pair[0]} vs {pair[1]}"
                        for pair in sample_set_overlap_results["identical elements, different order"]
                    ]
                )
            )
        if sample_set_overlap_results["partly overlapping"]:
            msg_lines.append(
                "Sample sets that are partly overlapping:\n"
                + "\n".join([f"{pair[0]} vs {pair[1]}" for pair in sample_set_overlap_results["partly overlapping"]])
            )
        full_msg = "\n\n".join(msg_lines)
        logger.error(full_msg)
        raise ValueError(full_msg)
    else:
        logger.info("No unsupported sample sets found. Proceeding with bcftools pipeline.")
        return


def _calculate_pairwise_overlap_types_for_sample_sets(sample_sets_dict: dict[tuple, list[str]]):
    """
    Analyze all pairwise overlap types between sample sets.
    Catches the case where two sets have the same elements but different order.

    Since dicts cannot have duplicate keys, there cannot be any completely overlapping sets in the input dict. Thus the else:continue condition should never occur.

    tuples are used in the input to maintain order within the sample sets, but tuples do not support '&' intersection operations. Thus they need to be converted with set()

    """
    keys = ["identical elements, different order", "partly overlapping", "non-overlapping"]
    sample_set_overlap_results = {key: [] for key in keys}
    sample_sets = list(sample_sets_dict.keys())
    logger.debug(f"Sample sets for overlap analysis: {sample_sets}")
    for set1, set2 in combinations(sample_sets, 2):
        set1_set = set(set1)
        set2_set = set(set2)
        if set1_set == set2_set and set1 != set2:
            overlap = "identical elements, different order"
        elif set1_set & set2_set:
            overlap = "partly overlapping"
        elif set1_set.isdisjoint(set2_set):
            overlap = "non-overlapping"
        else:
            continue
        sample_set_overlap_results[overlap].append((set1, set2))
    return sample_set_overlap_results


def _check_that_file_versions_match_dimensions_index(
    vcf_dimensions_data: dict,
    latest_versions_of_bucket_files: dict[str, str],
    metadata_result,
) -> None:
    """
    Ensure that the VCF dimensions index is up to date with the latest versions of the VCF files.
    """
    # Build lookup dict: filename -> s3_version_id
    vcf_lookup = {entry["vcf_file_s3_key"]: entry["s3_version_id"] for entry in vcf_dimensions_data["vcf_files"]}

    for file in metadata_result.unique_filenames:
        file_version_ID = latest_versions_of_bucket_files.get(file, "null")

        if file not in vcf_lookup or vcf_lookup[file] != file_version_ID:
            logger.error(f"VCF dimensions are outdated for file: {file}")
            raise ValueError(
                "The VCF dimensions file is not up to date with the VCF files in the project. "
                "Please run 'divbase-cli dimensions update --project <project_name>' and then submit the query again."
            )


def _record_task_metrics(task_metrics: TaskMetrics) -> None:
    """
    Helper function that records task metrics to the cache and updates Prometheus gauges
    so that they can be served by the worker metrics server.

    Note 1: There is a distinction between None (not captured) and 0.0 (ran but used no resources).
    None values are not stored in the cache or shown in Prometheus, while 0.0 values are.

    Note 2: task_cpu_seconds is the ARTIFICIAL SUM (python_overhead + bcftools),
    while task_python_overhead_cpu_seconds is the ACTUAL MEASURED value.
    """

    for field in dataclasses.fields(task_metrics):
        field_name = field.name
        if field_name in ("job_id", "task_name"):
            continue

        value = getattr(task_metrics, field_name)
        if value is not None:
            store_task_metric_in_cache(field_name, task_metrics.job_id, task_metrics.task_name, value)

    update_prometheus_gauges_from_cache()


# Import cron_tasks at the end to register all periodic tasks with the app. This avoids timing and circular import issues.
from divbase_api.worker import cron_tasks  # noqa: E402, F401
