"""
Frontend routes for user to manage their personal access tokens (PATs).

All routes here should rely on get_current_user_from_cookie dependency to ensure user is logged in.
"""

import logging
from datetime import datetime, timedelta, timezone
from zoneinfo import ZoneInfo

from fastapi import APIRouter, BackgroundTasks, Depends, Request, status
from fastapi.responses import HTMLResponse, RedirectResponse
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.personal_access_tokens import (
    create_personal_access_token,
    get_users_personal_access_tokens,
    soft_delete_personal_access_token,
)
from divbase_api.crud.projects import create_user_project_responses, get_user_projects_with_roles, has_required_role
from divbase_api.db import get_db
from divbase_api.deps import get_current_user_from_cookie
from divbase_api.exceptions import PATDuplicateNameError, PATLimitExceededError
from divbase_api.frontend_routes.core import templates
from divbase_api.models.projects import ProjectRoles
from divbase_api.models.users import UserDB
from divbase_api.services.email_sender import send_pat_created_email, send_pat_revoked_email
from divbase_lib.api_schemas.personal_access_tokens import PATPermissions

fr_pat_router = APIRouter()

logger = logging.getLogger(__name__)


@fr_pat_router.get("/", response_class=HTMLResponse, status_code=status.HTTP_200_OK)
async def list_pats_endpoint(
    request: Request,
    current_user: UserDB = Depends(get_current_user_from_cookie),
    db: AsyncSession = Depends(get_db),
    success: str | None = None,
    error: str | None = None,
):
    """Render the user's personal access tokens list page."""
    pats = await get_users_personal_access_tokens(db=db, user_id=current_user.id)
    return templates.TemplateResponse(
        request=request,
        name="pats_pages/personal_access_tokens.html",
        context={
            "current_user": current_user,
            "pats": pats,
            "now": datetime.now(timezone.utc),
            "success": success,
            "error": error,
        },
    )


@fr_pat_router.get("/new", response_class=HTMLResponse, status_code=status.HTTP_200_OK)
async def new_pat_form_endpoint(
    request: Request,
    current_user: UserDB = Depends(get_current_user_from_cookie),
    db: AsyncSession = Depends(get_db),
):
    """Render the new PAT creation form."""
    user_projects = await get_user_projects_with_roles(db=db, user_id=current_user.id)
    projects = create_user_project_responses(user_projects)
    return templates.TemplateResponse(
        request=request,
        name="pats_pages/new_personal_access_token.html",
        context={
            "current_user": current_user,
            "projects": projects,
        },
    )


@fr_pat_router.post("/new", response_class=HTMLResponse)
async def create_pat_endpoint(
    request: Request,
    background_tasks: BackgroundTasks,
    current_user: UserDB = Depends(get_current_user_from_cookie),
    db: AsyncSession = Depends(get_db),
):
    """
    Handle new PAT form submission.

    Unlike in other routes, the form data is parsed manually rather than as keyword args in fn signature.
    This is due to the complexity of the form data.
    """
    form = await request.form()
    name = str(form.get("name", "")).strip()
    description = str(form.get("description", "")).strip() or None
    expires_in_days = str(form.get("expires_in_days", "")).strip()
    scope_task_history = form.get("scope_task_history") == "on"
    scope_whoami = form.get("scope_whoami") == "on"
    project_access_mode = str(form.get("project_access_mode", "all"))

    user_projects = await get_user_projects_with_roles(db=db, user_id=current_user.id)
    projects = create_user_project_responses(user_projects)

    # Collect the user's per-project roles for re-rendering on error
    form_projects: dict[str, str] = {}
    for project in projects:
        role = str(form.get(f"project_{project.id}", "")).strip()
        if role:
            form_projects[str(project.id)] = role

    def form_error(message: str):
        return templates.TemplateResponse(
            request=request,
            name="pats_pages/new_personal_access_token.html",
            context={
                "error": message,
                "current_user": current_user,
                "projects": projects,
                "form_name": name,
                "form_description": description,
                "form_expires_in_days": expires_in_days,
                "form_scope_task_history": scope_task_history,
                "form_scope_whoami": scope_whoami,
                "form_project_access_specific": project_access_mode == "specific",
                "form_projects": form_projects,
            },
            status_code=status.HTTP_400_BAD_REQUEST,
        )

    if not name:
        return form_error("Token name is required.")
    if len(name) > 100:
        return form_error("Token name must be 100 characters or fewer.")
    if not expires_in_days:
        return form_error("Please select an expiration period.")

    expires_at_dt: datetime | None = None
    if expires_in_days != "never":
        expires_at_dt = datetime.now(timezone.utc) + timedelta(days=int(expires_in_days))

    if project_access_mode == "specific":
        user_role_by_project = {str(p.id): p.user_role for p in projects}
        for project_id_str, pat_role_str in form_projects.items():
            user_membership_role = user_role_by_project.get(project_id_str)
            if user_membership_role and not has_required_role(user_membership_role, ProjectRoles(pat_role_str)):
                project_name = next((p.name for p in projects if str(p.id) == project_id_str), project_id_str)
                return form_error(
                    f"Token role for '{project_name}' cannot exceed your current role ({user_membership_role})."
                )

    # pat_permissions=None means all user permissions.
    pat_permissions: PATPermissions | None = None
    if scope_task_history or scope_whoami or project_access_mode == "specific":
        pat_permissions = PATPermissions(
            task_history=scope_task_history,
            whoami=scope_whoami,
            all_projects=project_access_mode == "all",
            projects=form_projects if project_access_mode == "specific" else {},
        )

    try:
        pat, raw_token = await create_personal_access_token(
            db=db,
            user_id=current_user.id,
            name=name,
            description=description,
            permissions=pat_permissions.model_dump() if pat_permissions else None,
            expires_at=expires_at_dt,
        )
    except PATLimitExceededError as e:
        return form_error(e.message)
    except PATDuplicateNameError as e:
        return form_error(e.message)

    if expires_at_dt:
        dt = expires_at_dt.astimezone(ZoneInfo("Europe/Stockholm"))
        expires_at_cet = dt.strftime("%Y-%m-%d %H:%M:%S %Z")
    else:
        expires_at_cet = None

    background_tasks.add_task(
        send_pat_created_email, email_to=current_user.email, pat_name=name, expires_at_cet=expires_at_cet
    )
    return templates.TemplateResponse(
        request=request,
        name="pats_pages/personal_access_token_created.html",
        context={
            "current_user": current_user,
            "pat": pat,
            "expires_at_cet": expires_at_cet,
            "raw_token": raw_token.get_secret_value(),
        },
        status_code=status.HTTP_201_CREATED,
    )


@fr_pat_router.post("/{pat_id}/revoke", response_class=HTMLResponse)
async def revoke_pat_endpoint(
    pat_id: int,
    request: Request,
    background_tasks: BackgroundTasks,
    current_user: UserDB = Depends(get_current_user_from_cookie),
    db: AsyncSession = Depends(get_db),
):
    """Revoke (soft-delete) one of the current user's personal access tokens."""
    revoked_pat_name = await soft_delete_personal_access_token(db=db, pat_id=pat_id, user_id=current_user.id)
    if not revoked_pat_name:
        return RedirectResponse(
            url="/pats?error=Token+not+found+or+could+not+be+revoked",
            status_code=status.HTTP_302_FOUND,
        )

    background_tasks.add_task(send_pat_revoked_email, email_to=current_user.email, pat_name=revoked_pat_name)
    return RedirectResponse(
        url="/pats?success=Token+successfully+revoked",
        status_code=status.HTTP_302_FOUND,
    )
