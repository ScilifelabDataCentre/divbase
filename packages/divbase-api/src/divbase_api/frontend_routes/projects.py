"""
Frontend routes for user to view and (if manager) manage their projects.

For these routes you'll likely want to use one of the following dependencies:
- 'get_current_user_from_cookie':
    Just ensures the user is logged in
- 'get_project_member_from_cookie':
    Ensures the user is a member of a specific project and gives details about the project, user and their project role.
    (You need to validate their role is sufficient for what they're trying to do though).
"""

import logging

from fastapi import APIRouter, Depends, Form, HTTPException, Request, status
from fastapi.responses import HTMLResponse, RedirectResponse
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.projects import (
    add_project_member,
    create_user_project_responses,
    get_project_members,
    get_user_projects_with_roles,
    has_required_role,
    project_member_response_from_db,
    remove_project_member,
    update_project_member_role,
)
from divbase_api.crud.users import get_user_by_email
from divbase_api.db import get_db
from divbase_api.deps import get_current_user_from_cookie, get_project_member_from_cookie
from divbase_api.exceptions import AuthorizationError
from divbase_api.frontend_routes.core import templates
from divbase_api.models.projects import ProjectDB, ProjectRoles
from divbase_api.models.users import UserDB
from divbase_api.schemas.projects import UserProjectResponse

logger = logging.getLogger(__name__)

fr_projects_router = APIRouter()


@fr_projects_router.get("/", response_class=HTMLResponse)
async def get_user_projects_endpoint(
    request: Request,
    current_user: UserDB = Depends(get_current_user_from_cookie),
    db: AsyncSession = Depends(get_db),
):
    """Render the user's projects page showing all projects they're a member of."""
    user_projects = await get_user_projects_with_roles(db=db, user_id=current_user.id)
    projects = create_user_project_responses(user_projects)

    return templates.TemplateResponse(
        request=request,
        name="project_pages/index.html",
        context={
            "request": request,
            "current_user": current_user,
            "projects": projects,
        },
    )


@fr_projects_router.get("/{project_id}", response_class=HTMLResponse)
async def get_project_detail_endpoint(
    project_id: int,
    request: Request,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member_from_cookie),
    db: AsyncSession = Depends(get_db),
):
    """Render the project detail page."""
    project, current_user, role = project_and_user_and_role

    if not has_required_role(role, ProjectRoles.READ):
        raise AuthorizationError("You don't have permission to access this project.")

    project_response = UserProjectResponse(
        name=project.name,
        description=project.description,
        bucket_name=project.bucket_name,
        id=project.id,
        is_active=project.is_active,
        storage_quota_bytes=project.storage_quota_bytes,
        storage_used_bytes=project.storage_used_bytes,
        user_role=role,
    )

    members_db_model = await get_project_members(db=db, project_id=project_id)
    members = [project_member_response_from_db(membership) for membership in members_db_model]

    return templates.TemplateResponse(
        request=request,
        name="project_pages/detail.html",
        context={
            "request": request,
            "current_user": current_user,
            "project": project_response,
            "members": members,
        },
    )


@fr_projects_router.post("/{project_id}/members/add", response_class=HTMLResponse)
async def add_project_member_endpoint(
    project_id: int,
    user_email: str = Form(...),
    role: str = Form(...),
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member_from_cookie),
    db: AsyncSession = Depends(get_db),
):
    """Add a member to the project via frontend form."""
    project, current_user, user_role = project_and_user_and_role

    if not has_required_role(user_role, ProjectRoles.MANAGE):
        raise AuthorizationError("You don't have permission to manage this project.")

    user_to_add = await get_user_by_email(db=db, email=user_email.lower())
    if not user_to_add:
        # TODO - custom exception.
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"User with email '{user_email}' not found. Make sure they have an account.",
        )

    await add_project_member(db=db, project_id=project_id, user_id=user_to_add.id, role=ProjectRoles(role))
    logger.info(f"User '{user_to_add.email}' added to project {project.name} by user: {current_user.email}")
    return RedirectResponse(url=f"/projects/{project_id}", status_code=status.HTTP_303_SEE_OTHER)


@fr_projects_router.post("/{project_id}/members/{user_id}/role", response_class=HTMLResponse)
async def update_member_role_endpoint(
    project_id: int,
    user_id: int,
    role: str = Form(...),
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member_from_cookie),
    db: AsyncSession = Depends(get_db),
):
    """Update a member's role via frontend form."""
    project, current_user, user_role = project_and_user_and_role

    if not has_required_role(user_role, ProjectRoles.MANAGE):
        raise AuthorizationError("You don't have permission to manage this project.")

    _prevent_self_modification(user_id=user_id, current_user=current_user)

    _ = await update_project_member_role(db=db, project_id=project_id, user_id=user_id, new_role=ProjectRoles(role))
    logger.info(f"User '{user_id}' roles was updated in project {project.name} by user: {current_user.email}")
    return RedirectResponse(url=f"/projects/{project_id}", status_code=status.HTTP_303_SEE_OTHER)


@fr_projects_router.post("/{project_id}/members/{user_id}/remove", response_class=HTMLResponse)
async def remove_member_endpoint(
    project_id: int,
    user_id: int,
    project_and_user_and_role: tuple[ProjectDB, UserDB, ProjectRoles] = Depends(get_project_member_from_cookie),
    db: AsyncSession = Depends(get_db),
):
    """Remove a member from the project via frontend form."""
    project, current_user, user_role = project_and_user_and_role

    if not has_required_role(user_role, ProjectRoles.MANAGE):
        raise AuthorizationError("You don't have permission to manage this project.")

    _prevent_self_modification(user_id=user_id, current_user=current_user)

    await remove_project_member(db=db, project_id=project.id, user_id=user_id)
    logger.info(f"User '{user_id}' was removed from project {project.name} by user: {current_user.email}")
    return RedirectResponse(url=f"/projects/{project_id}", status_code=status.HTTP_303_SEE_OTHER)


def _prevent_self_modification(user_id: int, current_user: UserDB) -> None:
    """Prevent a user from modifying their own membership status in a project."""
    if user_id == current_user.id:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="You cannot remove or modify your own membership role in a project",
        )
