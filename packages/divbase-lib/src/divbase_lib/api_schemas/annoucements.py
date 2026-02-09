"""
Schemas for announcements.
"""

from pydantic import BaseModel, Field


class AnnouncementResponse(BaseModel):
    """Response model for an announcement returned by the API."""

    heading: str = Field(..., description="Title of the announcement.")
    message: str | None = Field(None, description="Detailed message of the announcement.")
    level: str = Field(
        ...,
        description="The announcement level, which can control styling of announcement. Possible values are: info, success, warning, danger.",
    )

    class Config:
        from_attributes = True
