"""
Schemas for working with S3 file operations.
"""

from pydantic import BaseModel


class DownloadOneObjectRequest(BaseModel):
    object_name: str
    version_id: str  # None means latest version


# Careful if trying to simplfy this logic, can a user download the same file at multiple versions?
class DownloadObjectsRequest(BaseModel):
    objects: list[DownloadOneObjectRequest]


class DownloadObjectResponse(BaseModel):
    """Response model for a single downloaded object with its pre-signed URL and version ID."""

    object_name: str
    pre_signed_url: str | None
    version_id: str | None  # None means latest version
