"""
Module responsible for taking pre-signed urls and using them to do file download and upload.

TODO: Consider adding retries, error handling, progress bars, etc.
"""

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import httpx
import stamina

from divbase_cli.user_auth import make_authenticated_request
from divbase_lib.api_schemas.s3 import (
    CompleteMultipartUploadRequest,
    CompleteMultipartUploadResponse,
    CreateMultipartUploadRequest,
    CreateMultipartUploadResponse,
    GetPresignedPartUrlsRequest,
    PreSignedDownloadResponse,
    PreSignedSinglePartUploadResponse,
    PresignedUploadPartUrlResponse,
    UploadedPart,
)
from divbase_lib.divbase_constants import (
    MAX_S3_API_BATCH_SIZE,
    S3_MULTIPART_CHUNK_SIZE,
)
from divbase_lib.exceptions import ChecksumVerificationError
from divbase_lib.s3_checksums import calculate_md5_checksum_for_chunk, verify_downloaded_checksum

logger = logging.getLogger(__name__)

# Used for multipart file transfers
MAX_CONCURRENCY = 10

MB = 1024 * 1024
# We can download in whatever chunk size we want,
# the checksum validation has to use the right chunk size though
# Uploads should use the same threshold and chunk size as divbase-api
# So those are defined in a shared lib.
MULTIPART_DOWNLOAD_THRESHOLD = 32 * MB
DOWNLOAD_CHUNK_SIZE = 8 * MB


@dataclass
class SuccessfulDownload:
    """Represents a successfully downloaded file."""

    file_path: Path
    object_name: str


@dataclass
class FailedDownload:
    """Represents a failed download attempt."""

    object_name: str
    file_path: Path
    exception: Exception


@dataclass
class DownloadOutcome:
    """Outcome of attempting to download multiple files."""

    successful: list[SuccessfulDownload]
    failed: list[FailedDownload]


def download_multiple_pre_signed_urls(
    pre_signed_urls: list[PreSignedDownloadResponse], verify_checksums: bool, download_dir: Path
) -> tuple[list[SuccessfulDownload], list[FailedDownload]]:
    """
    Download files using pre-signed URLs.
    Returns a tuple of both the successful and failed downloads.
    """
    successful_downloads, failed_downloads = [], []
    with httpx.Client() as client:
        for obj in pre_signed_urls:
            output_file_path = download_dir / obj.name
            object_name = obj.name
            try:
                result = _download_single_pre_signed_url(
                    httpx_client=client,
                    pre_signed_url=obj.pre_signed_url,
                    verify_checksums=verify_checksums,
                    output_file_path=output_file_path,
                    object_name=object_name,
                )
            except httpx.HTTPError as err:
                result = FailedDownload(object_name=object_name, file_path=output_file_path, exception=err)

            if isinstance(result, SuccessfulDownload):
                successful_downloads.append(result)
            else:
                failed_downloads.append(result)

    return successful_downloads, failed_downloads


def _download_single_pre_signed_url(
    httpx_client: httpx.Client, pre_signed_url: str, verify_checksums: bool, output_file_path: Path, object_name: str
) -> SuccessfulDownload | FailedDownload:
    """
    Download a single file using a pre-signed URL.
    If the file is large enough, we download in chunks.
    Helper function, do not call directly from outside this module.
    """
    content_length, server_checksum = _get_content_length_and_checksum(
        httpx_client=httpx_client, pre_signed_url=pre_signed_url
    )

    if content_length < MULTIPART_DOWNLOAD_THRESHOLD:
        _perform_singlepart_download(
            httpx_client=httpx_client,
            pre_signed_url=pre_signed_url,
            output_file_path=output_file_path,
        )

    else:
        logger.info(f"Starting multipart download for large file '{object_name}' of size {content_length} bytes.")
        _perform_multipart_download(httpx_client, pre_signed_url, output_file_path, content_length)

    if verify_checksums:
        try:
            verify_downloaded_checksum(file_path=output_file_path, expected_checksum=server_checksum)
        except ChecksumVerificationError as err:
            output_file_path.unlink()
            return FailedDownload(object_name=object_name, file_path=output_file_path, exception=err)

    return SuccessfulDownload(file_path=output_file_path, object_name=object_name)


def _get_content_length_and_checksum(httpx_client: httpx.Client, pre_signed_url: str) -> tuple[int, str]:
    """
    "HEAD" a pre-signed download URL to get it's content length and checksum.

    As you can't HEAD a presigned GET, we do a GET with a Range header to only get the first byte.
    Otherwise would have to be given a separate pre-signed HEAD url to do this.
    """
    with httpx_client.stream("GET", pre_signed_url, headers={"Range": "bytes=0-0"}) as head_response:
        head_response.raise_for_status()
        content_range = head_response.headers["Content-Range"]
        # format is "bytes 0-0/12345"
        content_length = int(content_range.split("/")[-1])
        server_checksum = head_response.headers["ETag"].strip('"')
    return content_length, server_checksum


def retry_only_on_retryable_http_errors(exc: Exception) -> bool:
    """
    Used by stamina's (library for retries) decorators to determine whether to retry the function or not.
    We avoid retrying on HTTPStatusError for 4xx errors as no point (e.g. 404 Not Found or 403 Forbidden etc...).
    """
    if isinstance(exc, httpx.HTTPStatusError):
        return exc.response.status_code >= 500

    # Want to retry on other HTTPError (parent of HTTPStatusError),
    # as this includes timeouts, connection errors, etc.
    return isinstance(exc, httpx.HTTPError)


@stamina.retry(on=retry_only_on_retryable_http_errors, attempts=3)
def _perform_singlepart_download(httpx_client: httpx.Client, pre_signed_url: str, output_file_path: Path) -> None:
    """Used on objects smaller than the multipart threshold cutoff"""
    with httpx_client.stream("GET", pre_signed_url) as response:
        response.raise_for_status()

        with open(output_file_path, "wb") as file:
            for chunk in response.iter_bytes(chunk_size=8192):
                file.write(chunk)


def _perform_multipart_download(httpx_client, pre_signed_url, output_file_path, content_length):
    """
    Download a large file in multiple chunks using range requests.

    As we write to the file concurrently, the file is first created with the correct size,
    and then each chunk is written to the correct position in the file.
    """
    with open(output_file_path, "wb") as f:
        f.seek(content_length - 1)
        f.write(b"\0")
    with ThreadPoolExecutor(max_workers=MAX_CONCURRENCY) as executor:
        futures = []
        for i in range(0, content_length, DOWNLOAD_CHUNK_SIZE):
            start = i
            end = min(i + DOWNLOAD_CHUNK_SIZE, content_length)
            futures.append(
                executor.submit(
                    _download_chunk,
                    httpx_client,
                    pre_signed_url,
                    start,
                    end,
                    output_file_path,
                )
            )

        for future in futures:
            future.result()


@stamina.retry(on=retry_only_on_retryable_http_errors, attempts=3)
def _download_chunk(client: httpx.Client, url: str, start: int, end: int, output_file_path: Path) -> None:
    """
    Downloads a specific range of bytes of a file (aka chunk),
    and writes it to the correct place in the file.
    """
    headers = {"Range": f"bytes={start}-{end - 1}"}
    with client.stream("GET", url, headers=headers) as response:
        response.raise_for_status()
        with open(output_file_path, "rb+") as f:
            f.seek(start)
            for chunk in response.iter_bytes():
                f.write(chunk)


@dataclass
class SuccessfulUpload:
    """Represents a successfully uploaded file."""

    file_path: Path
    object_name: str


@dataclass
class FailedUpload:
    """Represents a failed upload attempt."""

    object_name: str
    file_path: Path
    exception: Exception


@dataclass
class UploadOutcome:
    """Outcome of attempting to upload multiple files."""

    successful: list[SuccessfulUpload]
    failed: list[FailedUpload]


def upload_multiple_singlepart_pre_signed_urls(
    pre_signed_urls: list[PreSignedSinglePartUploadResponse], all_files: list[Path]
) -> UploadOutcome:
    """
    Upload singlepart files using pre-signed PUT URLs.
    Returns a UploadResults object containing the results of the upload attempts.
    """
    file_map = {file.name: file for file in all_files}

    successful_uploads, failed_uploads = [], []
    with httpx.Client() as client:
        for obj in pre_signed_urls:
            result = _upload_one_singlepart_pre_signed_url(
                httpx_client=client,
                pre_signed_url=obj.pre_signed_url,
                file_path=file_map[obj.name],
                object_name=obj.name,
                headers=obj.put_headers,
            )

            if isinstance(result, SuccessfulUpload):
                successful_uploads.append(result)
            else:
                failed_uploads.append(result)

    return UploadOutcome(successful=successful_uploads, failed=failed_uploads)


def _upload_one_singlepart_pre_signed_url(
    httpx_client: httpx.Client,
    pre_signed_url: str,
    file_path: Path,
    object_name: str,
    headers: dict[str, str],
) -> SuccessfulUpload | FailedUpload:
    """
    Upload one singlepart file to S3 using a pre-signed PUT URL.
    Helper function, do not call directly from outside this module.
    """
    with open(file_path, "rb") as file:
        try:
            response = httpx_client.put(pre_signed_url, content=file, headers=headers)
            response.raise_for_status()
        except httpx.HTTPError as err:
            return FailedUpload(object_name=object_name, file_path=file_path, exception=err)

    return SuccessfulUpload(file_path=file_path, object_name=object_name)


### multipart upload logic below ###


def perform_multipart_upload(
    project_name: str,
    divbase_base_url: str,
    file_path: Path,
    safe_mode: bool,
) -> SuccessfulUpload | FailedUpload:
    """
    Manages the entire multi-part upload process for a single file.
    Protocol as follows: TODO

    # TODO - error handling, which leads to abort the upload
    """
    object_name = file_path.name
    file_size = file_path.stat().st_size

    # 1. Create multipart upload
    create_request = CreateMultipartUploadRequest(
        name=object_name,
        content_length=file_size,
        part_size=S3_MULTIPART_CHUNK_SIZE,
    )
    response = make_authenticated_request(
        method="POST",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/s3/upload/multi-part/create?project_name={project_name}",
        json=create_request.model_dump(),
    )
    object_data = CreateMultipartUploadResponse(**response.json())

    # 2. Upload each part in batches as divbase server limits how many part urls it will give at once
    parts_to_request = list(range(1, object_data.number_of_parts + 1))

    uploaded_parts: list[UploadedPart] = []
    for i in range(0, len(parts_to_request), MAX_S3_API_BATCH_SIZE):
        part_batch_numbers = parts_to_request[i : i + MAX_S3_API_BATCH_SIZE]
        part_urls = _get_part_urls(
            project_name=project_name,
            divbase_base_url=divbase_base_url,
            object_name=object_name,
            upload_id=object_data.upload_id,
            part_numbers=part_batch_numbers,
            file_path=file_path,
            safe_mode=safe_mode,
        )
        batch_uploads = _upload_parts(part_urls=part_urls, file_path=file_path)
        uploaded_parts.extend(batch_uploads)

    # 3. Complete multipart upload
    complete_request_body = CompleteMultipartUploadRequest(
        name=object_name,
        upload_id=object_data.upload_id,
        parts=uploaded_parts,
    )
    response = make_authenticated_request(
        method="POST",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/s3/upload/multi-part/complete?project_name={project_name}",
        json=complete_request_body.model_dump(),
    )
    completed_upload = CompleteMultipartUploadResponse(**response.json())
    return SuccessfulUpload(file_path=file_path, object_name=completed_upload.name)

    # 4. TODO - implement error handling, abort upload and checksum validation for final object - is it needed?


def _get_part_urls(
    project_name: str,
    divbase_base_url: str,
    object_name: str,
    upload_id: str,
    part_numbers: list[int],
    file_path: Path,
    safe_mode: bool,
) -> list[PresignedUploadPartUrlResponse]:
    """
    Gets up to 100 pre-signed URLs (from divbase server) for uploading parts of a large file to S3.

    Not responsible for uploading the parts, just getting the URLs.
    """
    md5_checksums = None
    if safe_mode:
        md5_checksums = []
        for part_num in part_numbers:
            checksum = calculate_md5_checksum_for_chunk(
                file_path=file_path,
                start_byte=(part_num - 1) * S3_MULTIPART_CHUNK_SIZE,
                chunk_size=S3_MULTIPART_CHUNK_SIZE,
            )
            md5_checksums.append(checksum)

    request_body = GetPresignedPartUrlsRequest(
        name=object_name,
        upload_id=upload_id,
        parts_range_start=part_numbers[0],
        parts_range_end=part_numbers[-1],
        md5_checksums=md5_checksums,
    )
    response = make_authenticated_request(
        method="POST",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/s3/upload/multi-part/part-urls?project_name={project_name}",
        json=request_body.model_dump(),
    )
    return [PresignedUploadPartUrlResponse(**item) for item in response.json()]


def _upload_parts(part_urls: list[PresignedUploadPartUrlResponse], file_path: Path) -> list[UploadedPart]:
    """Uploads a batch of parts in parallel and returns their ETag info."""
    completed_parts = []
    with ThreadPoolExecutor(max_workers=MAX_CONCURRENCY) as executor:
        future_to_part = {executor.submit(_upload_chunk, part=part, file_path=file_path): part for part in part_urls}
        for future in as_completed(future_to_part):
            part_number, etag = future.result()
            completed_parts.append(UploadedPart(part_number=part_number, etag=etag))
    return completed_parts


@stamina.retry(on=retry_only_on_retryable_http_errors, attempts=3)
def _upload_chunk(part: PresignedUploadPartUrlResponse, file_path: Path) -> tuple[int, str]:
    """Uploads a single chunk of a file to a pre-signed URL and returns its part number and ETag."""

    start_byte = (part.part_number - 1) * S3_MULTIPART_CHUNK_SIZE
    with open(file_path, "rb") as f:
        f.seek(start_byte)
        data_to_upload = f.read(S3_MULTIPART_CHUNK_SIZE)

    with httpx.Client() as client:
        response = client.put(
            part.pre_signed_url,
            content=data_to_upload,
            headers=part.headers,
            timeout=httpx.Timeout(5.0, write=30.0),
        )
        response.raise_for_status()
        # ETag is returned with quotes, which must be stripped prior to comparison
        etag = response.headers["ETag"].strip('"')
        return part.part_number, etag
