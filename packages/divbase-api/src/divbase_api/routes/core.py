"""
Core API routes for divbase, including health checks and announcements.

Note that unlike every other API route these routes are not behind authentication...
"""

import logging

from fastapi import APIRouter, Depends, status
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.announcements import get_active_announcements
from divbase_api.db import get_db
from divbase_api.models.announcements import AnnouncementTarget
from divbase_lib.api_schemas.announcements import AnnouncementResponse

logger = logging.getLogger(__name__)

core_router = APIRouter()


@core_router.get("/health", status_code=status.HTTP_200_OK, response_model=dict[str, str])
def health():
    """Basic health check endpoint for the server."""
    return {"status": "ok"}


@core_router.get("/announcements", status_code=status.HTTP_200_OK, response_model=list[AnnouncementResponse])
async def announcements(db: AsyncSession = Depends(get_db)):
    """Returns active announcements for CLI users from the server."""
    return await get_active_announcements(db=db, target=AnnouncementTarget.CLI)
