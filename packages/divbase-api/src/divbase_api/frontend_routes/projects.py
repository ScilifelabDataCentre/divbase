"""
Frontend routes for user to view and (if manager) manage their projects.

These routes will return Template Responses.
"""

from fastapi import APIRouter, Depends, HTTPException, Request, status
from fastapi.responses import HTMLResponse
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.projects import (
    create_user_project_responses,
    get_project_with_user_role,
    get_user_projects_with_roles,
)
from divbase_api.db import get_db
from divbase_api.deps import get_current_user_from_cookie
from divbase_api.frontend_routes.core import templates
from divbase_api.models.users import UserDB
from divbase_api.schemas.projects import UserProjectResponse
from divbase_api.schemas.users import UserResponse

fr_projects_router = APIRouter()


@fr_projects_router.get("/", response_class=HTMLResponse)
async def user_projects_endpoint(
    request: Request,
    current_user: UserDB = Depends(get_current_user_from_cookie),
    db: AsyncSession = Depends(get_db),
):
    """Render the user's projects page showing all projects they're a member of."""
    user_projects = await get_user_projects_with_roles(db=db, user_id=current_user.id)
    user_project_responses = create_user_project_responses(user_projects)

    return templates.TemplateResponse(
        request=request,
        name="project_pages/index.html",
        context={
            "request": request,
            "current_user": UserResponse.model_validate(current_user),
            "projects": [UserProjectResponse.model_validate(project) for project in user_project_responses],
        },
    )


@fr_projects_router.get("/{project_id}", response_class=HTMLResponse)
async def project_detail_endpoint(
    project_id: int,
    request: Request,
    current_user: UserDB = Depends(get_current_user_from_cookie),
    db: AsyncSession = Depends(get_db),
):
    """Render the project detail page."""
    project, user_role = await get_project_with_user_role(db=db, project_id=project_id, user_id=current_user.id)

    if not project or not user_role:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Project not found or you don't have access")

    project_response = UserProjectResponse(
        name=project.name,
        description=project.description,
        bucket_name=project.bucket_name,
        id=project.id,
        is_active=project.is_active,
        storage_quota_bytes=project.storage_quota_bytes,
        storage_used_bytes=project.storage_used_bytes,
        user_role=user_role,
    )

    return templates.TemplateResponse(
        request=request,
        name="project_pages/detail.html",
        context={
            "request": request,
            "current_user": UserResponse.model_validate(current_user),
            "project": project_response,
        },
    )
