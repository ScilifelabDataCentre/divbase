import asyncio
import dataclasses
import logging
import os
import re
from itertools import combinations
from pathlib import Path

import httpx
from celery import Celery

from divbase_lib.exceptions import NoVCFFilesFoundError
from divbase_lib.queries import BCFToolsInput, BcftoolsQueryManager, run_sidecar_metadata_query
from divbase_lib.s3_client import S3FileManager, create_s3_file_manager
from divbase_lib.vcf_dimension_indexing import VCFDimensionCalculator

logger = logging.getLogger(__name__)

BROKER_URL = os.environ.get("CELERY_BROKER_URL", "pyamqp://guest@localhost//")
RESULT_BACKEND = os.environ.get("CELERY_RESULT_BACKEND", "redis://localhost:6379/0")

# Worker service account credentials
WORKER_SERVICE_EMAIL = os.environ.get("WORKER_SERVICE_EMAIL", "NOT_SET")
WORKER_SERVICE_PASSWORD = os.environ.get("WORKER_SERVICE_PASSWORD", "NOT_SET")
DIVBASE_API_URL = os.environ.get("DIVBASE_API_URL", "http://fastapi:8000")

app = Celery("divbase_worker", broker=BROKER_URL, backend=RESULT_BACKEND)

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


def _get_worker_access_token() -> str:
    """
    Get access token for worker service account.

    This follows the same pattern as the CLI login in user_auth.py.
    Token is fetched fresh for each task to avoid expiration issues.
    """
    if WORKER_SERVICE_EMAIL == "NOT_SET" or WORKER_SERVICE_PASSWORD == "NOT_SET":
        raise ValueError(
            "Worker service account credentials not set. "
            "Please set WORKER_SERVICE_EMAIL and WORKER_SERVICE_PASSWORD environment variables."
        )

    try:
        response = httpx.post(
            f"{DIVBASE_API_URL}/api/v1/auth/login",
            data={
                "grant_type": "password",
                "username": WORKER_SERVICE_EMAIL,
                "password": WORKER_SERVICE_PASSWORD,
            },
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            timeout=10.0,
        )
        response.raise_for_status()

        data = response.json()
        return data["access_token"]

    except httpx.HTTPError as e:
        logger.error(f"Failed to authenticate worker service account: {e}")
        raise


def _make_authenticated_api_request(method: str, endpoint: str, **kwargs) -> httpx.Response:
    """
    Make an authenticated request to the DivBase API.

    This follows the same pattern as CLI's make_authenticated_request in user_auth.py.

    Args:
        method: HTTP method (GET, POST, DELETE, etc.)
        endpoint: API endpoint path (e.g., "/api/v1/vcf-dimensions/projects/1")
        **kwargs: Additional arguments passed to httpx.request()

    Returns:
        httpx.Response: API response
    """
    access_token = _get_worker_access_token()

    headers = kwargs.pop("headers", {})
    headers["Authorization"] = f"Bearer {access_token}"

    url = f"{DIVBASE_API_URL}{endpoint}"

    try:
        response = httpx.request(method, url, headers=headers, timeout=30.0, **kwargs)
        response.raise_for_status()
        return response
    except httpx.HTTPError as e:
        logger.error(f"API request failed: {method} {endpoint} - {e}")
        raise


@app.task(name="tasks.sample_metadata_query", tags=["quick"])
def sample_metadata_query_task(tsv_filter: str, metadata_tsv_name: str, bucket_name: str) -> dict:
    """
    Run a sample metadata query task as a Celery task.
    """

    task_id = sample_metadata_query_task.request.id

    try:
        s3_file_manager = create_s3_file_manager(url="http://minio:9000")

        try:
            response = _make_authenticated_api_request(
                "GET", f"/api/v1/vcf-dimensions/lookup/project-by-bucket/{bucket_name}"
            )
            project_data = response.json()
            project_id = project_data["project_id"]
        except Exception as e:
            logger.error(f"Failed to get project info for bucket '{bucket_name}': {e}")
            return {"status": "error", "error": str(e), "type": "ProjectLookupError", "task_id": task_id}

        metadata_path = download_sample_metadata(
            metadata_tsv_name=metadata_tsv_name, bucket_name=bucket_name, s3_file_manager=s3_file_manager
        )

        try:
            response = _make_authenticated_api_request("GET", f"/api/v1/vcf-dimensions/list/project/{project_id}")
            vcf_dimensions_data = response.json()
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 404:
                logger.error(f"No VCF metadata found for project {project_id}")
                return {
                    "status": "error",
                    "error": f"No VCF dimensions indexed for project '{bucket_name}'. Please run 'divbase-cli dimensions update --project {bucket_name}' first.",
                    "type": "VCFDimensionsMissingError",
                    "task_id": task_id,
                }
            else:
                logger.error(f"Failed to get VCF dimensions: {e}")
                return {"status": "error", "error": str(e), "type": "VCFDimensionsAPIError", "task_id": task_id}

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

        # Return results (dataclass converted to dict for Celery serialization)
        result = dataclasses.asdict(metadata_result)
        result["status"] = "completed"
        result["task_id"] = task_id

        logger.info(
            f"Metadata query completed: {len(metadata_result.unique_sample_ids)} samples "
            f"mapped to {len(metadata_result.unique_filenames)} VCF files"
        )

        return result

    except Exception as e:
        logger.error(f"Error in sample_metadata_query_task: {str(e)}")
        return {
            "status": "error",
            "error": str(e),
            "type": type(e).__name__,
            "task_id": task_id,
        }


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

    try:
        response = _make_authenticated_api_request(
            "GET", f"/api/v1/vcf-dimensions/lookup/project-by-bucket/{bucket_name}"
        )
        project_data = response.json()
        project_id = project_data["project_id"]
    except Exception as e:
        logger.error(f"Failed to get project info for bucket '{bucket_name}': {e}")
        return {"status": "error", "error": str(e), "task_id": task_id}

    try:
        response = _make_authenticated_api_request("GET", f"/api/v1/vcf-dimensions/list/project/{project_id}")
        vcf_dimensions_data = response.json()

        if not vcf_dimensions_data.get("vcf_files"):
            error_msg = f"No VCF dimensions indexed for project '{bucket_name}'. Please run 'divbase-cli dimensions update --project {bucket_name}' first."
            logger.error(error_msg)
            return {"status": "error", "error": error_msg, "type": "VCFDimensionsMissingError", "task_id": task_id}

    except httpx.HTTPStatusError as e:
        if e.response.status_code == 404:
            error_msg = f"No VCF dimensions indexed for project '{bucket_name}'. Please run 'divbase-cli dimensions update --project {bucket_name}' first."
            logger.error(error_msg)
            return {"status": "error", "error": error_msg, "type": "VCFDimensionsMissingError", "task_id": task_id}
        else:
            logger.error(f"Failed to get VCF dimensions: {e}")
            return {"status": "error", "error": str(e), "type": "VCFDimensionsAPIError", "task_id": task_id}

    latest_versions_of_bucket_files = s3_file_manager.latest_version_of_all_files(bucket_name=bucket_name)

    metadata_path = download_sample_metadata(
        metadata_tsv_name=metadata_tsv_name, bucket_name=bucket_name, s3_file_manager=s3_file_manager
    )

    metadata_result = run_sidecar_metadata_query(
        file=metadata_path,
        filter_string=tsv_filter,
        project_id=project_id,
        vcf_dimensions_data=vcf_dimensions_data,
    )

    check_that_file_versions_match_dimensions_index(
        vcf_dimensions_data=vcf_dimensions_data,
        latest_versions_of_bucket_files=latest_versions_of_bucket_files,
        metadata_result=metadata_result,
    )

    files_to_download = metadata_result.unique_filenames
    sample_and_filename_subset = metadata_result.sample_and_filename_subset

    if "view -r" in command:
        files_to_download = check_for_unnecessary_files_for_region_query(
            command=command,
            files_to_download=metadata_result.unique_filenames,
            vcf_dimensions_data=vcf_dimensions_data,
        )

        sample_and_filename_subset = [
            entry for entry in metadata_result.sample_and_filename_subset if entry["Filename"] in files_to_download
        ]

    check_if_samples_can_be_combined_with_bcftools(
        files_to_download=files_to_download,
        vcf_dimensions_data=vcf_dimensions_data,
    )

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
        output_file = BcftoolsQueryManager().execute_pipe(command, bcftools_inputs, task_id)
    except Exception as e:
        logger.error(f"Error in bcftools task: {str(e)}")
        return {"status": "error", "error": str(e), "task_id": task_id}

    upload_results_file(output_file=Path(output_file), bucket_name=bucket_name, s3_file_manager=s3_file_manager)
    delete_job_files_from_worker(vcf_paths=files_to_download, metadata_path=metadata_path, output_file=output_file)

    return {"status": "completed", "output_file": output_file, "submitter": user_name}


@app.task(name="tasks.update_vcf_dimensions_task")
def update_vcf_dimensions_task(bucket_name: str, user_name: str = "Default User"):
    """
    Update VCF dimensions in the database for the specified bucket.

    This is a sync Celery task that wraps async operations via asyncio.run().
    """
    task_id = update_vcf_dimensions_task.request.id
    return asyncio.run(_update_vcf_dimensions_async(bucket_name, user_name, task_id))


async def _update_vcf_dimensions_async(bucket_name: str, user_name: str, task_id: str):
    """Async implementation - all database operations are async.

    Update the VCF dimensions file for the specified bucket. Ensures that the index only covers files in the bucket at the time of task execution.

    Re: vcfs_deleted_from_bucket_since_last_indexing: since both lists it is dependent on are already calculated,
    set difference is a faster operation than a new list comp.
    """
    s3_file_manager = create_s3_file_manager(url="http://minio:9000")

    all_files = s3_file_manager.list_files(bucket_name=bucket_name)
    vcf_files = [file for file in all_files if file.endswith(".vcf") or file.endswith(".vcf.gz")]

    if not vcf_files:
        raise NoVCFFilesFoundError(
            f"VCF dimensions file could not be generated since no VCF files were found in the project: {bucket_name}. Please upload at least one VCF file and run this command again."
        )

    try:
        response = _make_authenticated_api_request(
            "GET", f"/api/v1/vcf-dimensions/lookup/project-by-bucket/{bucket_name}"
        )
        project_data = response.json()
        project_id = project_data["project_id"]
    except Exception as e:
        logger.error(f"Failed to get project info for bucket '{bucket_name}': {e}")
        return {"status": "error", "error": str(e), "task_id": task_id}

    try:
        response = _make_authenticated_api_request("GET", f"/api/v1/vcf-dimensions/list/project/{project_id}")
        existing_metadata = response.json()
        already_indexed_vcfs = {
            entry["vcf_file_s3_key"]: entry["s3_version_id"] for entry in existing_metadata.get("vcf_files", [])
        }
    except httpx.HTTPStatusError as e:
        if e.response.status_code == 404:
            already_indexed_vcfs = {}
        else:
            logger.error(f"Failed to get existing VCF metadata: {e}")
            return {"status": "error", "error": str(e), "task_id": task_id}

    try:
        response = _make_authenticated_api_request("GET", f"/api/v1/vcf-dimensions/list-skipped/project/{project_id}")
        skipped_data = response.json()
        already_skipped_vcfs = {
            entry["vcf_file_s3_key"]: entry["s3_version_id"] for entry in skipped_data.get("skipped_files", [])
        }
    except httpx.HTTPStatusError as e:
        if e.response.status_code == 404:
            already_skipped_vcfs = {}
        else:
            logger.error(f"Failed to get skipped VCFs: {e}")
            already_skipped_vcfs = {}

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

    _ = download_vcf_files(
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

                response = _make_authenticated_api_request(
                    "POST", "/api/v1/vcf-dimensions/create-skipped", json=skipped_vcf_data
                )

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

            response = _make_authenticated_api_request("POST", "/api/v1/vcf-dimensions/create", json=vcf_metadata_data)

            files_indexed_by_this_job.append(file)
            logger.info(f"Indexed VCF metadata for: {file}")

        except Exception as e:
            logger.error(f"Error indexing {file}: {str(e)}")
            return {"status": "error", "error": str(e), "task_id": task_id}

    vcfs_deleted_from_bucket_since_last_indexing = list(set(already_indexed_vcfs) - set(vcf_files))
    skipped_deleted_from_bucket = list(set(already_skipped_vcfs) - set(vcf_files))

    if vcfs_deleted_from_bucket_since_last_indexing:
        for file in vcfs_deleted_from_bucket_since_last_indexing:
            try:
                response = _make_authenticated_api_request(
                    "DELETE", f"/api/v1/vcf-dimensions/delete/project/{project_id}/file/{file}"
                )
                logger.info(f"Deleted VCF metadata for removed file: {file}")
            except Exception as e:
                logger.error(f"Failed to delete VCF metadata for {file}: {e}")

    if skipped_deleted_from_bucket:
        for file in skipped_deleted_from_bucket:
            try:
                _make_authenticated_api_request(
                    "DELETE", f"/api/v1/vcf-dimensions/delete-skipped/project/{project_id}/file/{file}"
                )
                logger.info(f"Deleted skipped VCF entry for removed file: {file}")
            except Exception as e:
                logger.error(f"Failed to delete skipped VCF entry for {file}: {e}")

    delete_job_files_from_worker(vcf_paths=non_indexed_vcfs)

    if files_indexed_by_this_job == []:
        files_indexed_by_this_job = ["None: no new VCF files or file versions were detected in the project."]
    if divbase_results_files_skipped_by_this_job == []:
        divbase_results_files_skipped_by_this_job = ["None: no DivBase-generated results were detected in the project."]
    return {
        "status": "completed",
        "submitter": user_name,
        "VCF files that were added to dimensions index by this job": files_indexed_by_this_job,
        "VCF files skipped by this job (previous DivBase-generated result VCFs)": divbase_results_files_skipped_by_this_job,
        "VCF files that have been deleted from the project and thus have been dropped from the index": vcfs_deleted_from_bucket_since_last_indexing,
    }


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
    logger.debug(f"Starting download of {len(files_to_download)} VCF file(s) from bucket '{bucket_name}'")

    objects = {file_name: None for file_name in files_to_download}
    downloaded_files = s3_file_manager.download_files(
        objects=objects,
        download_dir=Path.cwd(),
        bucket_name=bucket_name,
    )

    logger.info(f"Downloaded VCF files: {[f.name for f in downloaded_files]}")

    return downloaded_files


def upload_results_file(output_file: Path, bucket_name: str, s3_file_manager: S3FileManager) -> None:
    """
    Upon completion of the task, upload the results file to the specified bucket.
    """
    _ = s3_file_manager.upload_files(
        to_upload={output_file.name: output_file},
        bucket_name=bucket_name,
    )


def check_for_unnecessary_files_for_region_query(
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


def delete_job_files_from_worker(
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


def check_if_samples_can_be_combined_with_bcftools(
    files_to_download: list[str],
    vcf_dimensions_data: dict,
) -> None:
    """
    Check if samples in VCF files can be combined with bcftools merge/concat.
    Raises ValueError if samples have incompatible overlaps.
    """
    if not vcf_dimensions_data.get("vcf_files"):
        raise ValueError("VCF dimensions data is missing or empty. Cannot check if samples can be combined.")

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

    sample_set_overlap_results = calculate_pairwise_overlap_types_for_sample_sets(sample_sets)
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


def calculate_pairwise_overlap_types_for_sample_sets(sample_sets_dict: dict[tuple, list[str]]):
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


def check_that_file_versions_match_dimensions_index(
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
