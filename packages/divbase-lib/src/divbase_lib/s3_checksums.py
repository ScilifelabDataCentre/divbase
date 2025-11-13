"""
Manages creation/validation of S3 object checksums for both file uploads and downloads.

We do not support multipart uploads/downloads at this time.
Single part uploads/downloads have a limit of 5GBs.
Docs: https://docs.netapp.com/us-en/storagegrid/s3/put-object.html
"""

import base64
import hashlib
import logging
from enum import StrEnum
from pathlib import Path
from typing import Iterator

from divbase_lib.exceptions import ChecksumVerificationError

logger = logging.getLogger(__name__)


def _read_file_chunks(file_path: Path, chunk_size: int = 8192) -> Iterator[bytes]:
    """Helper function to read a file in 'chunk_size' sized chunks."""

    with file_path.open(mode="rb") as infile:
        yield from iter(lambda: infile.read(chunk_size), b"")


class MD5CheckSumFormat(StrEnum):
    HEX = "hex"
    BASE64 = "base64"


def calculate_md5_checksum(file_path: Path, output_format: MD5CheckSumFormat) -> str:
    """
    Calculate the MD5 checksum of a file.
    Returns the checksum in either hex-encoded (lowercase) or base64-encoded format.
    """
    md5_hash = hashlib.md5()

    for chunk in _read_file_chunks(file_path):
        md5_hash.update(chunk)

    if output_format == MD5CheckSumFormat.HEX:
        return md5_hash.hexdigest()
    elif output_format == MD5CheckSumFormat.BASE64:
        return base64.b64encode(md5_hash.digest()).decode("utf-8")
    else:
        raise ValueError(f"Unknown output format: {output_format}")


def verify_downloaded_checksum(file_path: Path, expected_checksum: str) -> None:
    """
    Verify a downloaded file against S3's ETag (MD5 checksum in hex format).
    """
    calculated_md5 = calculate_md5_checksum(file_path=file_path, output_format=MD5CheckSumFormat.HEX)
    if calculated_md5 != expected_checksum:
        raise ChecksumVerificationError(expected_checksum=expected_checksum, calculated_checksum=calculated_md5)


def convert_checksum_hex_to_base64(hex_checksum: str) -> str:
    """
    Convert a hex-encoded MD5 checksum to base64-encoded format.
    """
    raw_bytes = bytes.fromhex(hex_checksum)
    base64_checksum = base64.b64encode(raw_bytes).decode("utf-8")
    return base64_checksum
