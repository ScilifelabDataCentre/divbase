"""
Frontend routes for core functionality.

E.g. home page, about page etc.

These routes will return Jinja2 Template Responses.
TODO - is it correct to say responsetype is html if we use Jinja2?
"""

from pathlib import Path

from fastapi import APIRouter, Depends, Request
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates

from divbase_api.deps import get_current_user_from_cookie_optional
from divbase_api.models.users import UserDB
from divbase_api.schemas.users import UserResponse

fr_core_router = APIRouter()

templates_dir = Path(__file__).parent.parent / "templates"
templates = Jinja2Templates(directory=templates_dir.resolve())


@fr_core_router.get("/", response_class=HTMLResponse)
async def get_home_page(request: Request, current_user: UserDB | None = Depends(get_current_user_from_cookie_optional)):
    """Render the home page."""
    if not current_user:
        return templates.TemplateResponse(
            request=request,
            name="index.html",
            context={"request": request, "user": None},
        )
    return templates.TemplateResponse(
        request=request,
        name="index.html",
        context={"request": request, "user": UserResponse.model_validate(current_user)},
    )
