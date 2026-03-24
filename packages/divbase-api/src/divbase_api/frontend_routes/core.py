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

from divbase_api import __version__ as divbase_version
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
templates.env.globals["divbase_version"] = divbase_version


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


def _simple_page(name: str, template: str):
    """Return a simple route handler that renders a basic static template for all core pages that don't require any extra context."""

    async def handler(
        request: Request,
        current_user: UserDB | None = Depends(get_current_user_from_cookie_optional),
    ):
        return templates.TemplateResponse(
            request=request,
            name=template,
            context={"current_user": current_user},
        )

    handler.__name__ = name
    return handler


fr_core_router.get("/about", response_class=HTMLResponse)(_simple_page("get_about", "core_pages/about.html"))
fr_core_router.get("/about/sv", response_class=HTMLResponse)(_simple_page("get_about_sv", "core_pages/about_sv.html"))
fr_core_router.get("/citation", response_class=HTMLResponse)(_simple_page("get_citation", "core_pages/citation.html"))
fr_core_router.get("/contact", response_class=HTMLResponse)(_simple_page("get_contact", "core_pages/contact.html"))
fr_core_router.get("/faqs", response_class=HTMLResponse)(_simple_page("get_faqs", "core_pages/faqs.html"))
fr_core_router.get("/terms", response_class=HTMLResponse)(_simple_page("get_terms", "core_pages/terms.html"))
fr_core_router.get("/privacy", response_class=HTMLResponse)(_simple_page("get_privacy", "core_pages/privacy.html"))
