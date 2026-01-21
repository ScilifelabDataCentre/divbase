"""
Schemas for DivBase's S3 API routes.

Pre-signed download URLs do not need to account for single vs multipart as this can be controlled by the client
using the HTTP range header when downloading (so you only need 1 pre-signed URL per object for download).

Pre-signed upload URLs need to account for single vs multipart uploads hence all the extra schemas below.
"""

from pydantic import BaseModel, Field

MB = 1024 * 1024


class DownloadObjectRequest(BaseModel):
    """Request model to upload a single object using a pre-signed URL."""

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
    part_size: int = Field(..., description="Size of each part in bytes", ge=8 * MB, le=64 * MB)


class CreateMultipartUploadResponse(BaseModel):
    """Response model to create a multipart upload using pre-signed URLs."""

    name: str = Field(..., description="Name of the object to be uploaded")
    upload_id: str = Field(..., description="Upload ID for the multipart upload")
    number_of_parts: int = Field(
        ..., description="Total number of parts required for the upload", ge=1, le=10000
    )  # TODO - could be determined client side?


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


# TODO - look into lifecycle rules for deleting incomplete multipart uploads.
class AbortMultipartUploadRequest(BaseModel):
    """Request model to abort a multipart upload and clean up parts."""

    name: str = Field(..., description="Name of the object being uploaded")
    upload_id: str = Field(..., description="Upload ID for the multipart upload to be aborted")


class AbortMultipartUploadResponse(BaseModel):
    """Response model to abort a multipart upload."""

    name: str = Field(..., description="Name of the object being uploaded")
    upload_id: str = Field(..., description="Upload ID for the multipart upload that was aborted")


class CheckFileExistsRequest(BaseModel):
    """Request model to check if a file already exists in the bucket (using the checksum)"""

    object_name: str
    md5_checksum: str


class ExistingFileResponse(BaseModel):
    """Response model for reporting a file that already exists in the bucket (using it's checksum)"""

    object_name: str
    md5_checksum: str
    matching_object_name: str | None
