"""
Unit tests for the creation of various S3 request objects.
Tests the validation logic correctly raises ValidationErrors when trying to build these objects.

NOTE: Tests validating object upload requests should validate both the
UploadSinglePartObjectRequest and CreateMultipartUploadRequest models.
"""

import pytest
from pydantic import ValidationError

from divbase_lib.api_schemas.s3 import (
    CreateMultipartUploadRequest,
    ListObjectsRequest,
    MakeDirectoriesRequest,
    UploadSinglePartObjectRequest,
)
from divbase_lib.divbase_constants import SUPPORTED_DIVBASE_FILE_TYPES, UNSUPPORTED_CHARACTERS_IN_FILENAMES

PATH_PREFIXES_WITH_DOT_SEGMENTS = ["../", "a/../", "a/./", "a/../../"]
PREFIX_ERROR_MATCH_STR = "cannot contain '.' or '..'"


@pytest.mark.parametrize("extension", SUPPORTED_DIVBASE_FILE_TYPES)
def test_upload_object_accepts_every_supported_file_type(extension):
    name = f"sample{extension}"
    single_part_request = UploadSinglePartObjectRequest(name=name, content_length=10, md5_hash=None)
    multi_part_request = CreateMultipartUploadRequest(name=name, content_length=10)
    assert name == single_part_request.name == multi_part_request.name


def test_upload_object_rejects_leading_or_trailing_slash():
    for name in ["/sample1.tsv", "sample1.tsv/"]:
        with pytest.raises(ValidationError, match="must not start or end with"):
            UploadSinglePartObjectRequest(name=name, content_length=10, md5_hash=None)
        with pytest.raises(ValidationError, match="must not start or end with"):
            CreateMultipartUploadRequest(name=name, content_length=10)


def test_upload_object_rejects_unsupported_file_extension():
    for extension in [".vcf", ".csv", ".exe"]:
        with pytest.raises(ValidationError, match="must be one of the supported file types"):
            UploadSinglePartObjectRequest(name=f"sample1{extension}", content_length=10, md5_hash=None)
        with pytest.raises(ValidationError, match="must be one of the supported file types"):
            CreateMultipartUploadRequest(name=f"sample1{extension}", content_length=10)


def test_upload_object_rejects_each_unsupported_character():
    for char in UNSUPPORTED_CHARACTERS_IN_FILENAMES:
        name = f"sample{char}.tsv"
        with pytest.raises(ValidationError, match="unsupported characters"):
            UploadSinglePartObjectRequest(name=name, content_length=10, md5_hash=None)
        with pytest.raises(ValidationError, match="unsupported characters"):
            CreateMultipartUploadRequest(name=name, content_length=10)


@pytest.mark.parametrize("path_prefix", PATH_PREFIXES_WITH_DOT_SEGMENTS)
def test_upload_object_rejects_double_dot_anywhere_in_the_path(path_prefix):
    file_name = f"{path_prefix}1a.tsv"

    with pytest.raises(ValidationError, match=PREFIX_ERROR_MATCH_STR):
        UploadSinglePartObjectRequest(name=file_name, content_length=10, md5_hash=None)
    with pytest.raises(ValidationError, match=PREFIX_ERROR_MATCH_STR):
        CreateMultipartUploadRequest(name=file_name, content_length=10)


def test_make_directories_request_accepts_valid_directories():
    """Test that MakeDirectoriesRequest accepts valid directory names."""
    directories = ["a/", "dir1/", "dir2/subdir/", "very/nested/sub/dir/"]
    request = MakeDirectoriesRequest(directories=directories)
    assert request.directories == directories

    upper_case_dirs = [d.upper() for d in directories]
    request = MakeDirectoriesRequest(directories=upper_case_dirs)
    assert request.directories == upper_case_dirs

    forward_slash_dirs = [d.rstrip("/") for d in directories]
    request = MakeDirectoriesRequest(directories=forward_slash_dirs)
    assert request.directories == directories  # forward slash gets normalised

    leading_slash_dirs = ["/" + d for d in directories]
    request = MakeDirectoriesRequest(directories=leading_slash_dirs)
    assert request.directories == directories  # leading slash gets stripped


def test_make_directories_request_normalizes_leading_and_trailing_slashes():
    possible_combos = ["/some/dir", "some/dir/", "/some/dir/", "some/dir"]
    desired_outcome = ["some/dir/"] * 4
    request = MakeDirectoriesRequest(directories=possible_combos)
    assert request.directories == desired_outcome


@pytest.mark.parametrize("bad_dir", ["", "/"])
def test_make_directories_request_rejects_empty_or_root_path(bad_dir):
    with pytest.raises(ValueError, match="cannot be empty or"):
        MakeDirectoriesRequest(directories=["good_dir/"] + [bad_dir])


def test_make_directories_request_rejects_each_unsupported_character():
    for char in UNSUPPORTED_CHARACTERS_IN_FILENAMES:
        with pytest.raises(ValidationError, match="unsupported characters"):
            MakeDirectoriesRequest(directories=[f"some{char}dir/"])


@pytest.mark.parametrize("path_prefix", PATH_PREFIXES_WITH_DOT_SEGMENTS)
def test_make_directories_request_rejects_double_dot_anywhere_in_the_path(path_prefix):
    with pytest.raises(ValidationError, match=PREFIX_ERROR_MATCH_STR):
        MakeDirectoriesRequest(directories=[path_prefix])


def test_list_objects_request_accepts_valid_prefixes():
    """Test that ListObjectsRequest accepts valid prefixes."""
    valid_prefixes = [None, "", "some/prefix", "some/dir/", "some/dir/file.txt"]
    for prefix in valid_prefixes:
        request = ListObjectsRequest(prefix=prefix, delimiter=None, next_token=None)
        if prefix == "":
            assert request.prefix is None
        else:
            assert request.prefix == prefix


def test_list_objects_request_rejects_each_unsupported_character():
    for char in UNSUPPORTED_CHARACTERS_IN_FILENAMES:
        with pytest.raises(ValidationError, match="unsupported characters"):
            ListObjectsRequest(prefix=f"some{char}prefix", delimiter=None, next_token=None)


@pytest.mark.parametrize("prefix", [".", "..", *PATH_PREFIXES_WITH_DOT_SEGMENTS])
def test_list_objects_request_rejects_double_dot_anywhere_in_the_path(prefix):
    with pytest.raises(ValidationError, match=PREFIX_ERROR_MATCH_STR):
        ListObjectsRequest(prefix=prefix, delimiter=None, next_token=None)
