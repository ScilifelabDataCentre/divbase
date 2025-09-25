"""
Pydantic schemas for admin-related operations.
"""

from pydantic import BaseModel, ConfigDict, Field


class SystemStatsResponse(BaseModel):
    """Response schema for system statistics. Used in admin dashboard."""

    total_users: int = Field(..., ge=0)
    admin_users: int = Field(..., ge=0)
    inactive_users: int = Field(..., ge=0)

    total_projects: int = Field(..., ge=0)
    active_projects: int = Field(..., ge=0)

    divbase_version: str = Field(..., description="DivBase version string")
    last_updated: str = Field(..., description="Last update")
    environment: str = Field(..., description="Deployment environment")

    model_config = ConfigDict(from_attributes=True)
