"""
Module responsible for taking pre-signed urls and using them to do file download and upload.

TODO: Consider adding retries, error handling, progress bars, etc.
"""

import logging
from dataclasses import dataclass
from pathlib import Path

import httpx

from divbase_cli.cli_exceptions import ChecksumVerificationError
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
    pre_signed_urls: list[dict], verify_checksums: bool, download_dir: Path
) -> DownloadOutcome:
    """
    Download files using pre-signed URLs.
    Returns a DownloadResults object containing all successful and failed downloads.
    """
    successful_downloads, failed_downloads = [], []
    with httpx.Client(timeout=30.0) as client:
        for obj in pre_signed_urls:
            object_name = obj["object_name"]
            pre_signed_url = obj["pre_signed_url"]
            out_file_path = download_dir / object_name

            with client.stream("GET", pre_signed_url) as response:
                try:
                    response.raise_for_status()
                except httpx.HTTPError as err:
                    failed_downloads.append(
                        FailedDownload(object_name=object_name, file_path=out_file_path, exception=err)
                    )
                    continue
                server_checksum = response.headers.get("ETag", "").strip('"')

                with open(out_file_path, "wb") as file:
                    for chunk in response.iter_bytes(chunk_size=8192):
                        file.write(chunk)

            if verify_checksums:
                try:
                    verify_downloaded_checksum(file_path=out_file_path, expected_checksum=server_checksum)
                    successful_downloads.append(SuccessfulDownload(file_path=out_file_path, object_name=object_name))

                except ChecksumVerificationError as err:
                    failed_downloads.append(
                        FailedDownload(object_name=object_name, file_path=out_file_path, exception=err)
                    )
            else:
                successful_downloads.append(SuccessfulDownload(file_path=out_file_path, object_name=object_name))

    return DownloadOutcome(successful=successful_downloads, failed=failed_downloads)


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


def upload_multiple_pre_signed_urls(pre_signed_urls: list[dict], all_files: list[Path]) -> UploadOutcome:
    """
    Upload files using pre-signed POST URLs.
    Returns a UploadResults object containing the results of the upload attempts.
    """
    file_map = {file.name: file for file in all_files}

    successful_uploads, failed_uploads = [], []

    with httpx.Client(timeout=30.0) as client:
        for obj in pre_signed_urls:
            object_name = obj["object_name"]
            post_url = obj["post_url"]
            fields = obj["fields"]

            file_path = file_map[object_name]

            with open(file_path, "rb") as file:
                files = {"file": (object_name, file, "application/octet-stream")}
                response = client.post(post_url, data=fields, files=files)

                try:
                    response.raise_for_status()
                    successful_uploads.append(SuccessfulUpload(file_path=file_path, object_name=object_name))
                except httpx.HTTPStatusError as err:
                    failed_uploads.append(FailedUpload(object_name=object_name, file_path=file_path, exception=err))

    return UploadOutcome(successful=successful_uploads, failed=failed_uploads)
