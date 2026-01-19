"""
Module responsible for taking pre-signed urls and using them to do file download and upload.

TODO: Consider adding retries, error handling, progress bars, etc.
"""

import logging
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import httpx

from divbase_lib.api_schemas.s3 import PreSignedDownloadResponse, PreSignedUploadResponse
from divbase_lib.exceptions import ChecksumVerificationError
from divbase_lib.s3_checksums import verify_downloaded_checksum

logger = logging.getLogger(__name__)

# Used for multipart file downloads
CHUNK_SIZE = 1024 * 1024 * 8  # 8MB
MAX_CONCURRENCY = 8
MULTIPART_DOWNLOAD_THRESHOLD = CHUNK_SIZE * 2


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
) -> DownloadOutcome:
    """
    Download files using pre-signed URLs.
    Returns a DownloadResults object containing all successful and failed downloads.
    """
    successful_downloads, failed_downloads = [], []
    with httpx.Client(timeout=30.0) as client:
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

    return DownloadOutcome(successful=successful_downloads, failed=failed_downloads)


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
        for i in range(0, content_length, CHUNK_SIZE):
            start = i
            end = min(i + CHUNK_SIZE, content_length)
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


def upload_multiple_pre_signed_urls(
    pre_signed_urls: list[PreSignedUploadResponse], all_files: list[Path]
) -> UploadOutcome:
    """
    Upload files using pre-signed PUT URLs.
    Returns a UploadResults object containing the results of the upload attempts.
    """
    file_map = {file.name: file for file in all_files}

    successful_uploads, failed_uploads = [], []
    with httpx.Client(timeout=30.0) as client:
        for obj in pre_signed_urls:
            result = _upload_single_pre_signed_url(
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


def _upload_single_pre_signed_url(
    httpx_client: httpx.Client,
    pre_signed_url: str,
    file_path: Path,
    object_name: str,
    headers: dict[str, str],
) -> SuccessfulUpload | FailedUpload:
    """
    Upload a single file using a pre-signed PUT URL.
    Helper function, do not call directly from outside this module.
    """
    with open(file_path, "rb") as file:
        try:
            response = httpx_client.put(pre_signed_url, content=file, headers=headers)
            response.raise_for_status()
        except httpx.HTTPError as err:
            return FailedUpload(object_name=object_name, file_path=file_path, exception=err)

    return SuccessfulUpload(file_path=file_path, object_name=object_name)
