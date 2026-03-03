"""
Crud operations for announcements

Announcements come primarily from the db model AnnouncementDB,
which stores announcements that can be displayed to users on the frontend and the cli.

Starlette admin will manage the creation/editing/deletion of announcements,
so this module only covers retrieving active announcements to be displayed on frontend or by CLI.
"""

from datetime import datetime, timezone

from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.api_config import settings
from divbase_api.models.announcements import AnnouncementDB, AnnouncementTarget
from divbase_api.services.validate_cli_versions import cli_update_available
from divbase_lib.api_schemas.announcements import AnnouncementResponse


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


def new_cli_version_announcement(user_cli_version: str) -> AnnouncementResponse | None:
    """
    Each request sent by CLI includes a header with user installed CLI version.
    We can use this to determine if we should add an announcement about a new version being available or not.
    """
    if not cli_update_available(cli_version=user_cli_version):
        return

    return AnnouncementResponse(
        heading="A new version of DivBase CLI is available",
        message=(
            f"You are using an outdated version of the DivBase CLI '{user_cli_version}'. "
            f"Please consider upgrading to the latest version '{settings.api.latest_cli_version}' for new features, bug fixes, and improved security.\n"
            "If you're not sure how to do that, you can find instructions on how to upgrade here: "
            f"{settings.api.mkdocs_site_url}/user-guides/installation"
        ),
        level="info",
    )
