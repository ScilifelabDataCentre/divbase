"""
Schemas for working with S3 file operations.
"""

from pydantic import BaseModel, Field


class DownloadOneObjectRequest(BaseModel):
    object_name: str = Field(..., description="Name of the object to be downloaded")
    version_id: str | None = Field(..., description="Version ID of the object, None if latest version")


# Careful if trying to simplfy this logic, can a user download the same file at multiple versions?
class DownloadObjectsRequest(BaseModel):
    objects: list[DownloadOneObjectRequest]


class PreSignedDownloadResponse(BaseModel):
    """Response model to download a single object using the pre-signed URL and (optionally) version ID."""

    object_name: str = Field(..., description="Name of the object to be downloaded")
    pre_signed_url: str = Field(..., description="Pre-signed URL for downloading the object")
    version_id: str | None = Field(..., description="Version ID of the object, None if latest version")


class PreSignedUploadResponse(BaseModel):
    """Response model to upload a single object using the pre-signed URL and field data."""

    object_name: str = Field(..., description="Name of the object to be uploaded")
    post_url: str = Field(..., description="Pre-signed URL to which the file should be uploaded")
    fields: dict = Field(..., description="Fields required for the POST request")
