"""
Pydantic schemas for project-related operations.
"""

from pydantic import BaseModel, ConfigDict, Field

from divbase_api.models.projects import ProjectRoles


class ProjectBase(BaseModel):
    name: str = Field(..., min_length=3, max_length=100)
    description: str | None = Field(None, min_length=1, max_length=1000)
    bucket_name: str = Field(..., min_length=3, max_length=63)

    model_config = ConfigDict(from_attributes=True, str_strip_whitespace=True)


class ProjectCreate(ProjectBase):
    storage_quota_bytes: int


class ProjectUpdate(BaseModel):
    name: str | None = Field(None, min_length=3, max_length=100)
    description: str | None = Field(None, min_length=1, max_length=1000)
    bucket_name: str = Field(..., min_length=3, max_length=63)
    storage_quota_bytes: int | None = None


class ProjectResponse(ProjectBase):
    """Response schema, aka returned by API endpoints."""

    id: int
    is_active: bool
    storage_quota_bytes: int  # TODO - more friendly format?
    storage_used_bytes: int


### ProjectMembership Schemas below ###


class ProjectMembershipCreate(BaseModel):
    role: ProjectRoles


class ProjectMembershipUpdate(BaseModel):
    role: ProjectRoles


class ProjectMembershipResponse(BaseModel):
    """Response schema, aka returned by API endpoints."""

    id: int
    user_id: int
    project_id: int
    role: ProjectRoles
