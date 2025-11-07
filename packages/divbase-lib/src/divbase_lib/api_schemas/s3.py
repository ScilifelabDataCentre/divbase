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
    md5_hash: str | None = Field(..., description="Optional MD5 hash of the object for integrity check")


class PreSignedUploadResponse(BaseModel):
    """Response model to upload a single object using the pre-signed URL and field data."""

    name: str = Field(..., description="Name of the object to be uploaded")
    post_url: str = Field(..., description="Pre-signed URL to which the file should be uploaded")
    fields: dict = Field(..., description="Fields required for the POST request")
