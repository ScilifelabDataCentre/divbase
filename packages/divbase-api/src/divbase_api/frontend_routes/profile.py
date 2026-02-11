"""
Frontend routes for user to manage their profile/account.

All routes here should rely on get_current_user_from_cookie dependency to ensure user is logged in.
"""

from fastapi import APIRouter, Depends, Form, Request, status
from fastapi.responses import HTMLResponse, RedirectResponse
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.projects import create_user_project_responses, get_user_projects_with_roles
from divbase_api.crud.users import update_user_profile
from divbase_api.db import get_db
from divbase_api.deps import get_current_user_from_cookie
from divbase_api.divbase_constants import SWEDISH_UNIVERSITIES
from divbase_api.frontend_routes.core import templates
from divbase_api.models.users import UserDB
from divbase_api.schemas.users import UserUpdate

fr_profile_router = APIRouter()


@fr_profile_router.get("/", response_class=HTMLResponse, status_code=status.HTTP_200_OK)
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


@fr_profile_router.get("/edit", response_class=HTMLResponse, status_code=status.HTTP_200_OK)
async def get_edit_user_profile_endpoint(
    request: Request,
    current_user: UserDB = Depends(get_current_user_from_cookie),
):
    """Render the edit user's profile page."""

    organisation_other = None
    organisation = current_user.organisation
    if current_user.organisation not in SWEDISH_UNIVERSITIES:
        organisation_other = current_user.organisation
        organisation = "Other"

    return templates.TemplateResponse(
        request=request,
        name="profile_pages/edit_profile.html",
        context={
            "request": request,
            "current_user": current_user,
            "swedish_universities": SWEDISH_UNIVERSITIES,
            # Pre-fill the form with the user's current information.
            # little awkward, but this logic is needed to handle the "Other" option for organisation.
            "organisation": organisation,
            "organisation_other": organisation_other,
            "organisation_role": current_user.organisation_role,
        },
    )


@fr_profile_router.post("/edit", response_class=HTMLResponse, status_code=status.HTTP_303_SEE_OTHER)
async def post_edit_user_profile_endpoint(
    request: Request,
    name: str = Form(...),
    current_user: UserDB = Depends(get_current_user_from_cookie),
    organisation: str = Form(...),
    organisation_other: str | None = Form(None),
    organisation_role: str = Form(...),
    db: AsyncSession = Depends(get_db),
):
    """
    Handle the submission of the edit user's profile form.

    In a similar manner to registration page, the organisation field is a dropdown with an "Other" option.
    If the user selects "Other", they must fill in the "organisation_other" field with their organisation name.
    We validate this here too.
    """
    if organisation == "Other":
        if not organisation_other:
            return templates.TemplateResponse(
                request=request,
                name="profile_pages/edit_profile.html",
                context={
                    "error": "Please specify your organisation.",
                    "request": request,
                    "current_user": current_user,
                    "swedish_universities": SWEDISH_UNIVERSITIES,
                    "name": name,
                    "organisation": organisation,
                    "organisation_other": organisation_other,
                    "organisation_role": organisation_role,
                },
            )
        else:
            organisation = organisation_other.strip()

    user_data = UserUpdate(name=name, organisation=organisation, organisation_role=organisation_role)
    _ = await update_user_profile(db=db, user_id=current_user.id, user_data=user_data)
    return RedirectResponse(url="/profile", status_code=status.HTTP_303_SEE_OTHER)
