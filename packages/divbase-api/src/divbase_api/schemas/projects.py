"""
Pydantic schemas for project-related operations.
"""

from pydantic import BaseModel, ConfigDict, Field

from divbase_api.models.projects import ProjectRoles
from divbase_lib.utils import format_file_size


class ProjectBase(BaseModel):
    name: str = Field(..., min_length=3, max_length=100)
    description: str | None = Field(None, min_length=0, max_length=1000)
    bucket_name: str = Field(..., min_length=3, max_length=63)

    model_config = ConfigDict(from_attributes=True, str_strip_whitespace=True)


class ProjectCreate(ProjectBase):
    storage_quota_bytes: int


class ProjectUpdate(BaseModel):
    name: str | None = Field(None, min_length=3, max_length=100)
    description: str | None = Field(None, min_length=0, max_length=1000)
    bucket_name: str | None = Field(None, min_length=3, max_length=63)
    storage_quota_bytes: int | None = None


class ProjectResponse(ProjectBase):
    """Response schema, aka returned by API endpoints."""

    id: int
    is_active: bool
    storage_quota_bytes: int
    storage_used_bytes: int

    @property
    def storage_used(self) -> str:
        """Convert storage used from bytes to a sensible value with units included."""
        return format_file_size(size_bytes=int(self.storage_used_bytes), decimals=1)

    @property
    def storage_quota(self) -> str:
        """Convert storage quota from bytes to a sensible value with units included."""
        return format_file_size(size_bytes=int(self.storage_quota_bytes), decimals=1)

    @property
    def storage_used_percent(self) -> int:
        """Calculate the percentage of storage used to nearest whole number."""
        if self.storage_quota_bytes == 0:
            return 0
        percent_used = round((self.storage_used_bytes / self.storage_quota_bytes) * 100)
        if percent_used == 0 and self.storage_used_bytes > 0:
            return 1  # avoid showing 0% when there is some storage used
        return percent_used

    model_config = ConfigDict(from_attributes=True, str_strip_whitespace=True)


class UserProjectResponse(ProjectResponse):
    """Response schema for projects that a user sees."""

    user_role: ProjectRoles


### ProjectMembership Schemas below ###


class ProjectMembershipCreate(BaseModel):
    role: ProjectRoles = Field(..., description="Role of the user in the project")


class ProjectMembershipUpdate(BaseModel):
    role: ProjectRoles = Field(..., description="Role of the user in the project")


class ProjectMembershipResponse(BaseModel):
    """Response schema, aka returned by API endpoints."""

    id: int
    user_id: int
    project_id: int
    role: ProjectRoles


class ProjectMemberResponse(BaseModel):
    """Response schema for project member with user details."""

    user_id: int
    user_name: str
    user_email: str
    user_is_active: bool
    role: ProjectRoles

    model_config = ConfigDict(from_attributes=True)
