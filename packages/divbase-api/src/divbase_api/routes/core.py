"""
Core API routes for divbase, including health checks and announcements.

Note that unlike every other API route these routes are not behind authentication...
"""

from fastapi import APIRouter, Depends, Request, status
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.announcements import get_active_announcements, new_cli_version_announcement
from divbase_api.db import get_db
from divbase_api.models.announcements import AnnouncementTarget
from divbase_lib.api_schemas.announcements import AnnouncementResponse
from divbase_lib.divbase_constants import CLI_VERSION_HEADER_KEY

core_router = APIRouter()


@core_router.get("/health", status_code=status.HTTP_200_OK, response_model=dict[str, str])
def health():
    """Basic health check endpoint for the server."""
    return {"status": "ok"}


@core_router.get("/announcements", status_code=status.HTTP_200_OK, response_model=list[AnnouncementResponse])
async def announcements(request: Request, db: AsyncSession = Depends(get_db)):
    """Returns active announcements for CLI users from the server."""
    db_announcements = await get_active_announcements(db=db, target=AnnouncementTarget.CLI)

    user_cli_version = request.headers.get(CLI_VERSION_HEADER_KEY)
    if not user_cli_version:
        return db_announcements

    new_version_announcement = new_cli_version_announcement(user_cli_version=user_cli_version)
    if new_version_announcement:
        return db_announcements + [new_version_announcement]
    return db_announcements
