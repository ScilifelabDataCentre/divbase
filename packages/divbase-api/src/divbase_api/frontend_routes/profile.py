"""
Frontend routes for user to manage their profile/account.

All routes here should rely on get_current_user_from_cookie dependency to ensure user is logged in.
"""

import logging

from fastapi import APIRouter, Depends, Form, Request, status
from fastapi.responses import HTMLResponse, RedirectResponse
from pydantic import ValidationError
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.api_constants import KNOWN_JOB_ROLES, SWEDISH_UNIVERSITIES
from divbase_api.crud.projects import create_user_project_responses, get_user_projects_with_roles
from divbase_api.crud.users import resolve_dropdown_form_input, update_user_profile
from divbase_api.db import get_db
from divbase_api.deps import get_current_user_from_cookie
from divbase_api.frontend_routes.core import templates
from divbase_api.models.users import UserDB
from divbase_api.schemas.users import UserUpdate

fr_profile_router = APIRouter()

logger = logging.getLogger(__name__)


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
    # little awkward, but this logic is needed to handle the "Other" option for
    # organisation and roles, and pre-fill each field appropriately on page load.
    organisation_other = None
    organisation = current_user.organisation
    if current_user.organisation not in SWEDISH_UNIVERSITIES:
        organisation_other = current_user.organisation
        organisation = "Other"

    role_other = None
    role = current_user.organisation_role
    if current_user.organisation_role not in KNOWN_JOB_ROLES:
        role_other = current_user.organisation_role
        role = "Other"

    return templates.TemplateResponse(
        request=request,
        name="profile_pages/edit_profile.html",
        context={
            "request": request,
            "current_user": current_user,
            "swedish_universities": SWEDISH_UNIVERSITIES,
            "known_job_roles": KNOWN_JOB_ROLES,
            "organisation": organisation,
            "organisation_other": organisation_other,
            "role": role,
            "role_other": role_other,
        },
    )


@fr_profile_router.post("/edit", response_class=HTMLResponse, status_code=status.HTTP_303_SEE_OTHER)
async def post_edit_user_profile_endpoint(
    request: Request,
    name: str = Form(...),
    current_user: UserDB = Depends(get_current_user_from_cookie),
    organisation: str = Form(...),
    organisation_other: str | None = Form(None),
    role: str = Form(...),
    role_other: str | None = Form(None),
    db: AsyncSession = Depends(get_db),
):
    """
    Handle the submission of the edit user's profile form.

    In a similar manner to registration page, the organisation field is a dropdown with an "Other" option.
    If the user selects "Other", they must fill in the "organisation_other" field with their organisation name.
    We validate this here too.
    """

    def edit_profile_failed_response(error_message: str):
        """Helper function to render the edit profile page with an error message."""
        return templates.TemplateResponse(
            request=request,
            name="profile_pages/edit_profile.html",
            context={
                "error": error_message,
                "request": request,
                "current_user": current_user,
                "swedish_universities": SWEDISH_UNIVERSITIES,
                "known_job_roles": KNOWN_JOB_ROLES,
                "name": name,
                "organisation": organisation,
                "organisation_other": organisation_other,
                "role": role,
                "role_other": role_other,
            },
        )

    resolved_organisation = resolve_dropdown_form_input(dropdown_value=organisation, other_value=organisation_other)
    if not resolved_organisation:
        return edit_profile_failed_response("Please specify your organisation, it must be at least 3 characters long.")

    resolved_role = resolve_dropdown_form_input(dropdown_value=role, other_value=role_other)
    if not resolved_role:
        return edit_profile_failed_response("Please specify your role, it must be at least 3 characters long.")

    try:
        user_data = UserUpdate(
            name=name.strip(),
            organisation=resolved_organisation.strip(),
            organisation_role=resolved_role.strip(),
        )
    except ValidationError as e:
        # This "should" not have happened, either:
        # Some mismatch in backend vs frontend validation logic or
        # someone bypassing client side validation (could be accidently or intentionally).
        logger.warning(f"User profile update failed backend validation for user_id: {current_user.id} - {e.errors()}")
        return edit_profile_failed_response(
            "Invalid input, please check your inputs match the required formats and try again."
        )

    _ = await update_user_profile(db=db, user_id=current_user.id, user_data=user_data)
    return RedirectResponse(url="/profile", status_code=status.HTTP_303_SEE_OTHER)
