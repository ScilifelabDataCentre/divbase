"""
Schemas for project versioning routes.

Project versions are the state of all files in a project's storage bucket at a given time point.
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
class ProjectBasicInfo(BaseModel):
    """Base model for describing a single project version, not for direct use in an endpoint."""

    name: str = Field(..., description="Version name")
    description: str | None = Field(..., description="Version description")


class AddVersionResponse(ProjectBasicInfo):
    """Response model for adding a version."""

    created_at: str = Field(..., description="ISO timestamp when version was created")


class ProjectVersionInfo(ProjectBasicInfo):
    """Basic information about a project version. You get a list of these when listing all versions in a project."""

    created_at: str = Field(..., description="ISO timestamp when version was created")
    is_deleted: bool = Field(..., description="Whether this version has been soft-deleted")


class ProjectVersionDetailResponse(ProjectBasicInfo):
    """Full information about a single project version, including the files at that version."""

    created_at: str = Field(..., description="ISO timestamp when version was created")
    is_deleted: bool = Field(..., description="Whether this version has been soft-deleted")
    files: dict[str, str] = Field(..., description="Mapping of file names to their version IDs")


class DeleteVersionResponse(BaseModel):
    """Response model for deleting a version."""

    name: str = Field(..., description="Name of the version that was deleted")
    already_deleted: bool = Field(
        False,
        description="Whether the version was already soft-deleted before this request",
    )
    date_deleted: str = Field(..., description="ISO timestamp of when the version was soft-deleted")
