"""
Schemas for DivBase's S3 API routes.

Pre-signed download URLs do not need to account for single vs multipart as this can be controlled by the client
using the HTTP range header when downloading (so you only need 1 pre-signed URL per object for download).

Pre-signed upload URLs need to account for single vs multipart uploads hence all the extra schemas below.
"""

from datetime import datetime

from pydantic import BaseModel, Field

from divbase_lib.divbase_constants import S3_MULTIPART_CHUNK_SIZE

MB = 1024 * 1024


## list objects models ##
class ListObjectsRequest(BaseModel):
    """Request model for listing objects in an S3 bucket."""

    prefix: str | None = Field(None, description="Optional prefix to filter objects by name.")
    next_token: str | None = Field(
        None, description="Token to continue listing files from the end of a previous request."
    )


class ObjectDetails(BaseModel):
    """Details about a single object in an S3 bucket."""

    name: str = Field(..., description="The name of the object in the bucket.")
    size: int = Field(..., description="The size of the object in bytes.")
    last_modified: datetime = Field(..., description="The date and time the object was last modified.")
    etag: str = Field(..., description="The ETag of the object, which is the MD5 checksum.")


class ListObjectsResponse(BaseModel):
    """Response model for listing objects in an S3 bucket."""

    objects: list[ObjectDetails] = Field(
        ..., description="A list of objects in the bucket.", min_length=0, max_length=1000
    )
    next_token: str | None = Field(
        None, description="Token for fetching the next page of results. If None, no more results."
    )


## file info models ##
class ObjectVersionInfo(BaseModel):
    """Detailed information about a single version of an S3 object."""

    version_id: str = Field(..., description="The version ID of the object.")
    last_modified: datetime = Field(..., description="The date and time the object version was last modified.")
    size: int = Field(..., description="The size of the object in bytes.")
    etag: str = Field(..., description="The ETag of the object, which is the MD5 checksum.")
    is_latest: bool = Field(..., description="Indicates if this is the latest version of the object.")


class ObjectInfoResponse(BaseModel):
    """Response model for detailed information about all versions of a single object stored in S3."""

    object_name: str = Field(..., description="The name of the object.")
    is_currently_deleted: bool = Field(..., description="True if the latest version of the object is a delete marker.")
    versions: list[ObjectVersionInfo] = Field(..., description="A list of all versions of the object.")


## download models ##
class DownloadObjectRequest(BaseModel):
    """Request model to download a single object using a pre-signed URL."""

    name: str = Field(..., description="Name of the object to be downloaded")
    version_id: str | None = Field(..., description="Version ID of the object, None if latest version")


class PreSignedDownloadResponse(BaseModel):
    """Response model to download a single object using the pre-signed URL and (optionally) version ID."""

    name: str = Field(..., description="Name of the object to be downloaded")
    pre_signed_url: str = Field(..., description="Pre-signed URL for downloading the object")
    version_id: str | None = Field(..., description="Version ID of the object, None if latest version")


### Single-part upload models ###
class UploadSinglePartObjectRequest(BaseModel):
    """Request model to upload a single object as a single part using a pre-signed URL."""

    name: str = Field(..., description="Name of the object to be uploaded")
    content_length: int = Field(..., description="Size of the file in bytes")
    md5_hash: str | None = Field(None, description="Optional MD5 hash of the object for integrity check")


class PreSignedSinglePartUploadResponse(BaseModel):
    """Response model to upload a single object as a single part using the pre-signed URL using PUT."""

    name: str = Field(..., description="Name of the object to be uploaded")
    pre_signed_url: str = Field(..., description="Pre-signed URL to which the file should be uploaded")
    put_headers: dict[str, str] = Field(..., description="Headers to be included in the PUT request")


### Multipart upload models ###
class CreateMultipartUploadRequest(BaseModel):
    """Request model to create a multipart upload using pre-signed URLs."""

    name: str = Field(..., description="Name of the object to be uploaded")
    content_length: int = Field(..., description="Size of the file in bytes")


class CreateMultipartUploadResponse(BaseModel):
    """Response model to create a multipart upload using pre-signed URLs."""

    name: str = Field(..., description="Name of the object to be uploaded")
    upload_id: str = Field(..., description="Upload ID for the multipart upload")
    number_of_parts: int = Field(..., description="Total number of parts required for the upload", ge=1, le=10000)
    part_size: int = Field(
        S3_MULTIPART_CHUNK_SIZE, description="Size of each part in bytes (the last part may be smaller)."
    )


class GetPresignedPartUrlsRequest(BaseModel):
    """
    Request model to get pre-signed URLs for multiple parts of a presigned multipart upload.

    You can request up to 100 parts at a time.
    Part number indexing is 1-based (with max allowed range: 1 to 10000).
    """

    name: str = Field(..., description="Name of the object to be uploaded")
    upload_id: str = Field(..., description="Upload ID for the multipart upload")
    parts_range_start: int = Field(..., description="Starting part number", ge=1, le=10000)
    parts_range_end: int = Field(..., description="Ending part number", ge=1, le=10000)
    md5_checksums: list[str] | None = Field(
        None, description="Optional list of MD5 checksums for each part to be uploaded"
    )


class PresignedUploadPartUrlResponse(BaseModel):
    """Response model for a pre-signed URL for a single part of a multipart upload."""

    part_number: int = Field(..., description="Part number", ge=1, le=10000)
    pre_signed_url: str = Field(..., description="Pre-signed URL for uploading this part")
    headers: dict[str, str] = Field(..., description="Headers to be included in the PUT request for this part")


class UploadedPart(BaseModel):
    """Model representing a part of an object that has been uploaded via multi-part upload."""

    part_number: int = Field(..., description="Part number", ge=1, le=10000)
    etag: str = Field(description="ETag returned by S3 after uploading the part")


class CompleteMultipartUploadRequest(BaseModel):
    """Request model to complete a multipart upload using pre-signed URLs."""

    name: str = Field(..., description="Name of the object to be uploaded")
    upload_id: str = Field(..., description="Upload ID for the multipart upload")
    parts: list[UploadedPart] = Field(..., description="List of parts that have been uploaded")


class CompleteMultipartUploadResponse(BaseModel):
    """Response model to complete a multipart upload using pre-signed URLs."""

    name: str = Field(..., description="Name of the object that was uploaded")
    version_id: str = Field(..., description="Version ID of the uploaded object")
    md5_hash: str = Field(..., description="MD5 hash of the uploaded object")


class AbortMultipartUploadRequest(BaseModel):
    """Request model to abort a multipart upload and clean up parts."""

    name: str = Field(..., description="Name of the object being uploaded")
    upload_id: str = Field(..., description="Upload ID for the multipart upload to be aborted")


class AbortMultipartUploadResponse(BaseModel):
    """Response model to abort a multipart upload."""

    name: str = Field(..., description="Name of the object being uploaded")
    upload_id: str = Field(..., description="Upload ID for the multipart upload that was aborted")


class RestoreObjectsResponse(BaseModel):
    """Response model for restoring soft-deleted objects in a bucket."""

    restored: list[str] = Field(
        ...,
        description="List of object names that were successfully restored, this includes objects that were already live",
    )
    not_restored: list[str] = Field(
        ...,
        description=(
            "List of object names that could not be processed.\n"
            "This could be due to several reasons:\n"
            "1. The object does not exist in the bucket (e.g., a typo in the name).\n"
            "2. The object was hard-deleted and is unrecoverable.\n"
            "3. An unexpected server error occurred during the restore attempt."
        ),
    )


## checksum models ##
class FileChecksumResponse(BaseModel):
    """Response model for reporting a file's checksum in the bucket."""

    object_name: str = Field(..., description="Name of the object in the bucket")
    md5_checksum: str = Field(..., description="MD5 checksum of the object in the bucket")
