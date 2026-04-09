"""
Frontend routes for user to manage their personal access tokens (PATs).

All routes here should rely on get_current_user_from_cookie dependency to ensure user is logged in.
"""

import logging
from datetime import datetime, timezone

from fastapi import APIRouter, Depends, Request, status
from fastapi.responses import HTMLResponse
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.personal_access_tokens import (
    get_users_personal_access_tokens,
)
from divbase_api.db import get_db
from divbase_api.deps import get_current_user_from_cookie
from divbase_api.frontend_routes.core import templates
from divbase_api.models.users import UserDB

fr_pat_router = APIRouter()

logger = logging.getLogger(__name__)


@fr_pat_router.get("/", response_class=HTMLResponse, status_code=status.HTTP_200_OK)
async def list_pats_endpoint(
    request: Request,
    current_user: UserDB = Depends(get_current_user_from_cookie),
    db: AsyncSession = Depends(get_db),
    success: str | None = None,
    error: str | None = None,
):
    """Render the user's personal access tokens list page."""
    pats = await get_users_personal_access_tokens(db=db, user_id=current_user.id)
    return templates.TemplateResponse(
        request=request,
        name="pats_pages/personal_access_tokens.html",
        context={
            "current_user": current_user,
            "pats": pats,
            "now": datetime.now(timezone.utc),
            "success": success,
            "error": error,
        },
    )
