"""
Frontend routes for Admin to manage DivBase projects

All routes here should rely on get_current_admin_user_from_cookie dependency to ensure user is admin.
"""

import logging

from fastapi import APIRouter, Depends, Request
from fastapi.responses import HTMLResponse
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.projects import get_all_projects
from divbase_api.db import get_db
from divbase_api.deps import get_current_admin_user_from_cookie
from divbase_api.frontend_routes.core import templates
from divbase_api.models.users import UserDB
from divbase_api.schemas.users import UserResponse

logger = logging.getLogger(__name__)

fr_admin_projects_router = APIRouter()


@fr_admin_projects_router.get("/", response_class=HTMLResponse)
async def admin_projects_endpoint(
    request: Request,
    current_admin: UserDB = Depends(get_current_admin_user_from_cookie),
    db: AsyncSession = Depends(get_db),
):
    """Render the admin's projects management page."""
    all_projects = await get_all_projects(db=db, limit=100)
    return templates.TemplateResponse(
        request=request,
        name="admin_pages/projects.html",
        context={
            "request": request,
            "current_user": UserResponse.model_validate(current_admin),
            "all_projects": all_projects,
        },
    )
