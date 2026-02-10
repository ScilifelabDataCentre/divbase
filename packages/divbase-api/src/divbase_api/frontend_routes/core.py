"""
Frontend routes for core pages (e.g. home page, about page etc...)

For these routes you will likely want to use the 'get_current_user_from_cookie_optional' dependency
to get the current user if they are logged in, but not require it.
"""

from pathlib import Path

from fastapi import APIRouter, Depends, Request
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.api_config import settings
from divbase_api.crud.announcements import get_active_announcements
from divbase_api.db import get_db
from divbase_api.deps import get_current_user_from_cookie_optional
from divbase_api.models.announcements import AnnouncementTarget
from divbase_api.models.users import UserDB

fr_core_router = APIRouter()

templates_dir = Path(__file__).parent.parent / "templates"
templates = Jinja2Templates(directory=templates_dir.resolve())
templates.env.globals["mkdocs_site_url"] = settings.api.mkdocs_site_url


@fr_core_router.get("/", response_class=HTMLResponse)
async def get_home_page(
    request: Request,
    db: AsyncSession = Depends(get_db),
    current_user: UserDB | None = Depends(get_current_user_from_cookie_optional),
):
    """Render the home page."""
    announcements = await get_active_announcements(db=db, target=AnnouncementTarget.WEB)
    return templates.TemplateResponse(
        request=request,
        name="index.html",
        context={"request": request, "current_user": current_user, "announcements": announcements},
    )
