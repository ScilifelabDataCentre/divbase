"""
Schemas for bucket versioning routes.
"""

from pydantic import BaseModel, Field


# Request Models
class CreateVersioningFileRequest(BaseModel):
    name: str = Field(..., description="Initial version name")
    description: str = Field(..., description="Initial version description")


class AddVersionRequest(BaseModel):
    name: str = Field(..., description="Name of the new version to add")
    description: str = Field("", description="Description of the new version")


class DeleteVersionRequest(BaseModel):
    version_name: str = Field(..., description="Name of the version to delete")


# Response Models


class CreateVersioningFileResponse(BaseModel):
    """Response model for creating versioning file."""

    name: str = Field(..., description="Name of the initial version created")
    description: str = Field(..., description="Description of the initial version created")


class AddVersionResponse(BaseModel):
    """Response model for adding a version."""

    name: str = Field(..., description="Name of the added version")
    description: str = Field(..., description="Description of the added version")


class BucketVersionDetail(BaseModel):
    """Full information about a single bucket version."""

    timestamp: str = Field(..., description="ISO timestamp when version was created")
    description: str = Field(..., description="Version description")
    files: dict[str, str] = Field(..., description="Mapping of file names to their version IDs")


class VersionListResponse(BaseModel):
    """Response model for listing all versions."""

    versions: dict[str, BucketVersionDetail] = Field(..., description="All versions in the bucket")


class FilesAtVersionResponse(BaseModel):
    """Response model for files at a specific version."""

    files: dict[str, str] = Field(..., description="Mapping of file names to their version IDs at a specific version")


class DeleteVersionResponse(BaseModel):
    deleted_version: str = Field(..., description="Name of the version that was deleted")
