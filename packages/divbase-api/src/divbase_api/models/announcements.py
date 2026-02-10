"""
Announcements DB Model.

Stores announcements that can be displayed to users on the frontend (banner) or via the CLI.
"""

from datetime import datetime
from enum import StrEnum

from sqlalchemy import DateTime, Enum, String
from sqlalchemy.orm import Mapped, mapped_column

from divbase_api.models.base import BaseDBModel


class AnnouncementTarget(StrEnum):
    """Possible targets for announcements."""

    BOTH = "both"
    CLI = "cli"
    WEB = "web"


class AnnouncementLevel(StrEnum):
    """
    Possible levels for announcements.
    These match bootstrap alert levels, so this will control the announcement styling on the frontend.
    """

    INFO = "info"
    SUCCESS = "success"
    WARNING = "warning"
    DANGER = "danger"


class AnnouncementDB(BaseDBModel):
    """
    DB Model for an announcement.

    id, created_at and updated_at are inherited from BaseDBModel.
    """

    __tablename__ = "announcement"

    heading: Mapped[str] = mapped_column(String(200), index=True)
    message: Mapped[str | None] = mapped_column(String(1000), nullable=True)
    target: Mapped[AnnouncementTarget] = mapped_column(Enum(AnnouncementTarget), index=True)
    level: Mapped[AnnouncementLevel] = mapped_column(Enum(AnnouncementLevel), index=True)
    auto_expire_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)

    def __repr__(self) -> str:
        return f"<AnnouncementDB id={self.id}, heading={self.heading}, target={self.target}, level={self.level}>"
