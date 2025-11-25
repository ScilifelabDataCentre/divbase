"""
Module responsible for taking pre-signed urls and using them to do file download and upload.

TODO: Consider adding retries, error handling, progress bars, etc.
"""

import logging
from dataclasses import dataclass
from pathlib import Path

import httpx

from divbase_lib.api_schemas.s3 import PreSignedDownloadResponse, PreSignedUploadResponse
from divbase_lib.exceptions import ChecksumVerificationError
from divbase_lib.s3_checksums import verify_downloaded_checksum

logger = logging.getLogger(__name__)


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
            result = _download_single_pre_signed_url(
                httpx_client=client,
                pre_signed_url=obj.pre_signed_url,
                verify_checksums=verify_checksums,
                output_file_path=download_dir / obj.name,
                object_name=obj.name,
            )
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
    Helper function, do not call directly from outside this module.
    """
    with httpx_client.stream("GET", pre_signed_url) as response:
        try:
            response.raise_for_status()
        except httpx.HTTPError as err:
            return FailedDownload(object_name=object_name, file_path=output_file_path, exception=err)

        server_checksum = response.headers.get("ETag", "").strip('"')

        with open(output_file_path, "wb") as file:
            for chunk in response.iter_bytes(chunk_size=8192):
                file.write(chunk)

    if verify_checksums:
        try:
            verify_downloaded_checksum(file_path=output_file_path, expected_checksum=server_checksum)
        except ChecksumVerificationError as err:
            return FailedDownload(object_name=object_name, file_path=output_file_path, exception=err)

    return SuccessfulDownload(file_path=output_file_path, object_name=object_name)


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
