"""
Frontend routes for user to manage their profile/account.

These routes will return Template Responses.

TODO: Currently only handle GET requests.
"""

from fastapi import APIRouter, Depends, Request
from fastapi.responses import HTMLResponse

from divbase_api.deps import get_current_user_from_cookie
from divbase_api.frontend_routes.core import templates
from divbase_api.models.users import UserDB
from divbase_api.schemas.users import UserResponse

fr_profile_router = APIRouter()


@fr_profile_router.get("/", response_class=HTMLResponse)
async def user_profile_endpoint(
    request: Request,
    current_user: UserDB = Depends(get_current_user_from_cookie),
):
    """Render the user's profile page with their information."""
    return templates.TemplateResponse(
        request=request,
        name="profile_pages/index.html",
        context={
            "request": request,
            "current_user": UserResponse.model_validate(current_user),
        },
    )


@fr_profile_router.get("/edit", response_class=HTMLResponse)
async def get_edit_user_profile_endpoint(
    request: Request,
    current_user: UserDB = Depends(get_current_user_from_cookie),
):
    """Render the edit user's profile page."""
    return templates.TemplateResponse(
        request=request,
        name="profile_pages/edit_profile.html",
        context={
            "request": request,
            "current_user": UserResponse.model_validate(current_user),
        },
    )
