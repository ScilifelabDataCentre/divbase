"""
Frontend routes for Admin to manage the site.

These routes will return Template Responses.

All routes here should rely on get_current_admin_user_from_cookie dependency to ensure user is admin.
"""

import logging

from fastapi import APIRouter, Depends, Request
from fastapi.responses import HTMLResponse
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.admin import get_system_stats
from divbase_api.db import get_db
from divbase_api.deps import get_current_admin_user_from_cookie
from divbase_api.frontend_routes.core import templates
from divbase_api.models.users import UserDB
from divbase_api.schemas.admin import SystemStatsResponse
from divbase_api.schemas.users import UserResponse

logger = logging.getLogger(__name__)

fr_admin_router = APIRouter()


@fr_admin_router.get("/", response_class=HTMLResponse)
async def admin_index_endpoint(
    request: Request,
    current_admin_user: UserDB = Depends(get_current_admin_user_from_cookie),
    db: AsyncSession = Depends(get_db),
):
    """Render the admin index page."""
    stats = await get_system_stats(db=db)
    return templates.TemplateResponse(
        request=request,
        name="admin_pages/index.html",
        context={
            "request": request,
            "user": UserResponse.model_validate(current_admin_user),
            "stats": SystemStatsResponse.model_validate(stats),
        },
    )
