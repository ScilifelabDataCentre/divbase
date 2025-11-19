"""
Schemas for DivBase's S3 API routes.
"""

from pydantic import BaseModel, Field


class DownloadObjectRequest(BaseModel):
    """Request model to upload a single object using a pre-signed URL."""

    name: str = Field(..., description="Name of the object to be downloaded")
    version_id: str | None = Field(..., description="Version ID of the object, None if latest version")


class PreSignedDownloadResponse(BaseModel):
    """Response model to download a single object using the pre-signed URL and (optionally) version ID."""

    name: str = Field(..., description="Name of the object to be downloaded")
    pre_signed_url: str = Field(..., description="Pre-signed URL for downloading the object")
    version_id: str | None = Field(..., description="Version ID of the object, None if latest version")


class UploadObjectRequest(BaseModel):
    """Request model to upload a single object using a pre-signed URL."""

    name: str = Field(..., description="Name of the object to be uploaded")
    content_length: int = Field(..., description="Size of the file in bytes")
    md5_hash: str | None = Field(None, description="Optional MD5 hash of the object for integrity check")


class PreSignedUploadResponse(BaseModel):
    """Response model to upload a single object using the pre-signed URL using PUT."""

    name: str = Field(..., description="Name of the object to be uploaded")
    pre_signed_url: str = Field(..., description="Pre-signed URL to which the file should be uploaded")
    put_headers: dict[str, str] = Field(..., description="Headers to be included in the PUT request")


class CheckFileExistsRequest(BaseModel):
    """Request model to check if a file already exists in the bucket (using the checksum)"""

    object_name: str
    md5_checksum: str


class ExistingFileResponse(BaseModel):
    """Response model for reporting a file that already exists in the bucket (using it's checksum)"""

    object_name: str
    md5_checksum: str
    matching_object_name: str | None
