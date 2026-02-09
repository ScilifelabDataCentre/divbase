"""
Crud operations on the accouncements table,
which stores announcements that can be displayed to users on the frontend and the cli.

Starlette admin will manage the creation/editing/deletion of announcements,
so this module only covers retrieving active announcements to be displayed on frontend or by CLI.
"""

from datetime import datetime, timezone

from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.models.annoucements import AnnouncementDB, AnnouncementTarget
from divbase_lib.api_schemas.annoucements import AnnouncementResponse


async def get_active_announcements(db: AsyncSession, target: AnnouncementTarget) -> list[AnnouncementResponse]:
    """Get active announcements for a given target (cli, web or both)."""

    if target == AnnouncementTarget.BOTH:
        raise ValueError("Target cannot be both when retrieving announcements. Please specify either cli or web.")

    stmt = (
        select(AnnouncementDB)
        .where((AnnouncementDB.target == target) | (AnnouncementDB.target == AnnouncementTarget.BOTH))
        .where((AnnouncementDB.auto_expire_at.is_(None)) | (AnnouncementDB.auto_expire_at > datetime.now(timezone.utc)))
    )
    result = await db.execute(stmt)
    announcements = result.scalars().all()
    return [AnnouncementResponse.model_validate(announcement) for announcement in announcements]
