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

# Chunk size used for multipart uploads to S3.
# When uploading large files, they are split into many chunks
# To validate the integrity of the downloaded file, you have to know the chunk size used.
S3_MULTIPART_CHUNK_SIZE = 32 * 1024 * 1024  # 32 MiB


class MD5CheckSumFormat(StrEnum):
    HEX = "hex"
    BASE64 = "base64"


def verify_downloaded_checksum(
    file_path: Path,
    expected_checksum: str,
) -> None:
    """
    Verify a downloaded file against it's S3's ETag.

    For files uploaded as single part, this is just the MD5 checksum in hex format.
    For files uploaded as multipart, this is a composite checksum of all the parts
    """
    if "-" in expected_checksum:
        calculated_checksum = calculate_composite_md5_s3_etag(file_path)
    else:
        calculated_checksum = calculate_md5_checksum(file_path=file_path, output_format=MD5CheckSumFormat.HEX)

    if calculated_checksum != expected_checksum:
        raise ChecksumVerificationError(expected_checksum=expected_checksum, calculated_checksum=calculated_checksum)


def _read_file_chunks(file_path: Path, chunk_size: int) -> Iterator[bytes]:
    """Helper function to read a file in 'chunk_size' sized chunks."""

    with file_path.open(mode="rb") as infile:
        yield from iter(lambda: infile.read(chunk_size), b"")


def calculate_md5_checksum(
    file_path: Path, output_format: MD5CheckSumFormat, chunk_size: int = S3_MULTIPART_CHUNK_SIZE
) -> str:
    """
    Calculate the MD5 checksum of a file.
    Returns the checksum in either hex-encoded (lowercase) or base64-encoded format.

    Used for:
    - Generating the "Content-MD5" header for S3 uploads (base64-encoded)
    - Verifying downloaded files against S3 ETag (hex-encoded)
        (only works for files uploaded as single part - not composite/multipart)
    """
    md5_hash = hashlib.md5()

    for chunk in _read_file_chunks(file_path=file_path, chunk_size=chunk_size):
        md5_hash.update(chunk)

    if output_format == MD5CheckSumFormat.HEX:
        return md5_hash.hexdigest()
    elif output_format == MD5CheckSumFormat.BASE64:
        return base64.b64encode(md5_hash.digest()).decode("utf-8")
    else:
        raise ValueError(f"Unknown output format: {output_format}")


def calculate_md5_checksum_for_chunk(file_path: Path, start_byte: int, chunk_size: int) -> str:
    """
    Calculate the base64-encoded MD5 checksum for a specific chunk of a file.
    S3 uses this checksum (Content-MD5 header) when uploading parts of a file.
    """
    md5_hash = hashlib.md5()
    with file_path.open("rb") as f:
        f.seek(start_byte)
        chunk = f.read(chunk_size)
        md5_hash.update(chunk)
    return base64.b64encode(md5_hash.digest()).decode("utf-8")


def calculate_composite_md5_s3_etag(
    file_path: Path,
    chunk_size: int = S3_MULTIPART_CHUNK_SIZE,
) -> str:
    """
    Calculate the composite ETag for a file that was uploaded via multipart upload to S3.
    This is used to validate the downloaded file's integrity.

    The process involves calculating the MD5 hash of each part, then combining these hashes to form a final ETag.
    So the part size used here must match the part size used during upload.
    """
    md5_digests = []
    part_count = 0

    for chunk in _read_file_chunks(file_path=file_path, chunk_size=chunk_size):
        md5_digests.append(hashlib.md5(chunk).digest())
        part_count += 1

    composite_hash = hashlib.md5(b"".join(md5_digests))
    return f"{composite_hash.hexdigest()}-{part_count}"


def convert_checksum_hex_to_base64(hex_checksum: str) -> str:
    """
    Convert a hex-encoded MD5 checksum to base64-encoded format.
    """
    raw_bytes = bytes.fromhex(hex_checksum)
    base64_checksum = base64.b64encode(raw_bytes).decode("utf-8")
    return base64_checksum
