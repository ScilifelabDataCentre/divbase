"""
Schemas for announcements.
"""

from pydantic import BaseModel


class AnnouncementResponse(BaseModel):
    """Response model for an announcement."""

    heading: str
    message: str | None
    level: str

    class Config:
        from_attributes = True
