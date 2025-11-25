"""
Frontend routes for user to manage their profile/account.

All routes here should rely on get_current_user_from_cookie dependency to ensure user is logged in.

TODO: (some of this relies on having email sending set up)
- Change password
- Delete account
- Change email
"""

from fastapi import APIRouter, Depends, Form, Request, status
from fastapi.responses import HTMLResponse, RedirectResponse
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.projects import create_user_project_responses, get_user_projects_with_roles
from divbase_api.crud.users import update_user_profile
from divbase_api.db import get_db
from divbase_api.deps import get_current_user_from_cookie
from divbase_api.frontend_routes.core import templates
from divbase_api.models.users import UserDB
from divbase_api.schemas.users import UserUpdate

fr_profile_router = APIRouter()


@fr_profile_router.get("/", response_class=HTMLResponse)
async def user_profile_endpoint(
    request: Request,
    current_user: UserDB = Depends(get_current_user_from_cookie),
    db: AsyncSession = Depends(get_db),
):
    """Render the user's profile page with their information."""
    user_projects = await get_user_projects_with_roles(db=db, user_id=current_user.id)
    projects = create_user_project_responses(user_projects)

    return templates.TemplateResponse(
        request=request,
        name="profile_pages/index.html",
        context={
            "request": request,
            "current_user": current_user,
            "projects": projects,
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
            "current_user": current_user,
        },
    )


@fr_profile_router.post("/edit", response_class=HTMLResponse)
async def post_edit_user_profile_endpoint(
    name: str = Form(...),
    current_user: UserDB = Depends(get_current_user_from_cookie),
    db: AsyncSession = Depends(get_db),
):
    """Handle the submission of the edit user's profile form."""
    user_data = UserUpdate(name=name)
    _ = await update_user_profile(db=db, user_id=current_user.id, user_data=user_data)
    return RedirectResponse(url="/profile", status_code=status.HTTP_303_SEE_OTHER)
