"""
Frontend routes for core pages (e.g. home page, about page etc...)

For these routes you will likely want to use the 'get_current_user_from_cookie_optional' dependency
to get the current user if they are logged in, but not require it.
"""

from pathlib import Path

from fastapi import APIRouter, Depends, Request
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates

from divbase_api.deps import get_current_user_from_cookie_optional
from divbase_api.models.users import UserDB

fr_core_router = APIRouter()

templates_dir = Path(__file__).parent.parent / "templates"
templates = Jinja2Templates(directory=templates_dir.resolve())


@fr_core_router.get("/", response_class=HTMLResponse)
async def get_home_page(request: Request, current_user: UserDB | None = Depends(get_current_user_from_cookie_optional)):
    """Render the home page."""
    return templates.TemplateResponse(
        request=request,
        name="index.html",
        context={"request": request, "current_user": current_user},
    )
