"""
Unit tests for the s3_checksums module.
"""

import base64
import hashlib
import random
from pathlib import Path

import pytest

from divbase_lib.divbase_constants import (
    S3_MULTIPART_UPLOAD_THRESHOLD,
)
from divbase_lib.exceptions import ChecksumVerificationError
from divbase_lib.s3_checksums import (
    MD5CheckSumFormat,
    calculate_composite_md5_s3_etag,
    calculate_md5_checksum,
    calculate_md5_checksum_for_chunk,
    convert_checksum_hex_to_base64,
    verify_downloaded_checksum,
)


@pytest.fixture
def small_test_file(tmp_path: Path) -> Path:
    """Creates a small file with known content."""
    file_path = tmp_path / "small_file.txt"
    file_path.write_text("This is a test file for checksums.")
    return file_path


@pytest.fixture
def large_test_file(tmp_path: Path) -> Path:
    """
    Creates a file larger than one multipart upload threshold with random but deterministic content.
    (So we can validate the checksum calculations)
    """
    file_path = tmp_path / "large_file.bin"
    one_mib = 1024 * 1024
    file_size = S3_MULTIPART_UPLOAD_THRESHOLD + one_mib

    random.seed(42)
    with open(file_path, "wb") as f:
        for _ in range(file_size // one_mib):
            f.write(random.randbytes(one_mib))

    assert file_path.stat().st_size == file_size
    return file_path


def test_calculate_md5_checksum(small_test_file: Path):
    """
    Tests MD5 checksum calculation in both hex and base64 formats.
    """
    # pre-determined checksum for "This is a test file for checksums."
    known_hex = "c5fdf96f527b6276e70d8aab9988eb24"
    known_base64 = "xf35b1J7YnbnDYqrmYjrJA=="

    hex_checksum = calculate_md5_checksum(small_test_file, MD5CheckSumFormat.HEX)
    base64_checksum = calculate_md5_checksum(small_test_file, MD5CheckSumFormat.BASE64)

    assert hex_checksum == known_hex
    assert base64_checksum == known_base64


def test_calculate_md5_checksum_for_chunk(small_test_file: Path):
    """
    Tests calculating a base64 MD5 checksum for a specific chunk of a file.
    """
    # (These parts have to match the content of the 'small_test_file' fixture)
    chunk1_content = b"This is a "
    chunk2_content = b"test file "
    base64_for_chunk2 = hashlib.md5(chunk2_content).digest()
    expected_checksum = base64.b64encode(base64_for_chunk2).decode("utf-8")

    calculated_checksum = calculate_md5_checksum_for_chunk(
        file_path=small_test_file, start_byte=len(chunk1_content), chunk_size=len(chunk2_content)
    )
    assert calculated_checksum == expected_checksum


def test_calculate_composite_md5_s3_etag(large_test_file: Path):
    """
    Tests the calculation of a composite ETag for a multipart file.
    Hard to test with known values since the file is random, so we just check the format.

    """
    known_etag = "733fce785e4a6c35e24c9a6ccf0e9d19-4"
    calculated_etag = calculate_composite_md5_s3_etag(large_test_file)
    assert calculated_etag == known_etag


def test_verify_downloaded_checksum_single_part_success(small_test_file: Path):
    """
    Tests successful verification for a single-part ETag.
    """
    correct_checksum = "c5fdf96f527b6276e70d8aab9988eb24"
    verify_downloaded_checksum(file_path=small_test_file, expected_checksum=correct_checksum)


def test_verify_downloaded_checksum_single_part_failure(small_test_file: Path):
    """
    Tests that ChecksumVerificationError is raised for a single-part mismatch.
    """
    incorrect_checksum = "incorrectchecksum123"
    with pytest.raises(ChecksumVerificationError):
        verify_downloaded_checksum(file_path=small_test_file, expected_checksum=incorrect_checksum)


def test_verify_downloaded_checksum_multi_part_success(large_test_file: Path):
    """
    Tests successful verification for a multi-part ETag.
    """
    correct_checksum = calculate_composite_md5_s3_etag(large_test_file)
    verify_downloaded_checksum(file_path=large_test_file, expected_checksum=correct_checksum)


def test_verify_downloaded_checksum_multi_part_failure(large_test_file: Path):
    """
    Tests that ChecksumVerificationError is raised for a multi-part mismatch.
    """
    incorrect_checksum = "incorrectcomposite-etag-2"
    with pytest.raises(ChecksumVerificationError):
        verify_downloaded_checksum(file_path=large_test_file, expected_checksum=incorrect_checksum)


def test_convert_checksum_hex_to_base64():
    """
    Tests the conversion from a hex checksum to its base64 representation.
    """
    hex_checksum = "a59336128d5045d419b1554978a8a378"
    expected_base64 = "pZM2Eo1QRdQZsVVJeKijeA=="
    converted = convert_checksum_hex_to_base64(hex_checksum)
    assert converted == expected_base64
