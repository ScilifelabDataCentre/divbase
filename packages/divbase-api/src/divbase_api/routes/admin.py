"""
Admin routes for basic tasks like creating users and projects.
This is only used in local development and testing to automatically setup test/dev data.

These routes are turned off in other enviroments.
We manage deployed instances of DivBase via the admin-panel (which uses starlette-admin).

All routes in here should depend on get_current_admin_user.
"""

import logging

from fastapi import APIRouter, BackgroundTasks, Depends, Query, status
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.projects import add_project_member, create_project
from divbase_api.crud.users import create_user
from divbase_api.db import get_db
from divbase_api.deps import get_current_admin_user
from divbase_api.models.projects import ProjectRoles
from divbase_api.models.users import UserDB
from divbase_api.schemas.projects import ProjectCreate, ProjectMembershipResponse, ProjectResponse
from divbase_api.schemas.users import UserCreate, UserResponse
from divbase_api.services.email_sender import send_test_email

logger = logging.getLogger(__name__)

admin_router = APIRouter()


@admin_router.post("/users/", response_model=UserResponse, status_code=status.HTTP_201_CREATED)
async def create_user_endpoint(
    user_data: UserCreate,
    is_admin: bool = Query(False, description="Set to true to create an admin user"),
    email_verified: bool = Query(False, description="Set to true to skip the email verification process"),
    db: AsyncSession = Depends(get_db),
    current_admin: UserDB = Depends(get_current_admin_user),
):
    """Create a new regular or admin user."""
    new_user = await create_user(db=db, user_data=user_data, is_admin=is_admin, email_verified=email_verified)
    logger.info(f"Admin user: {current_admin.email} created a new user: {new_user.email}")
    return new_user


@admin_router.post("/projects", response_model=ProjectResponse, status_code=status.HTTP_201_CREATED)
async def create_project_endpoint(
    proj_data: ProjectCreate,
    db: AsyncSession = Depends(get_db),
    current_admin: UserDB = Depends(get_current_admin_user),
):
    new_project = await create_project(db=db, proj_data=proj_data)
    logger.info(f"Admin user: {current_admin.email} created a new project: {new_project.name}")
    return new_project


@admin_router.post(
    "/projects/{project_id}/members/{user_id}",
    response_model=ProjectMembershipResponse,
    status_code=status.HTTP_201_CREATED,
)
async def add_project_member_endpoint(
    project_id: int,
    user_id: int,
    role: ProjectRoles,
    db: AsyncSession = Depends(get_db),
    current_admin: UserDB = Depends(get_current_admin_user),
):
    membership = await add_project_member(db=db, project_id=project_id, user_id=user_id, role=role)
    logger.info(
        f"Admin user: {current_admin.email} added user with id: {user_id} to project with id: {project_id} with role: {role}"
    )

    return ProjectMembershipResponse(
        id=membership.id,
        project_id=membership.project_id,
        user_id=membership.user_id,
        role=membership.role,
    )


@admin_router.post("/test-email/{email_to}", status_code=status.HTTP_200_OK)
async def test_email_endpoint(
    email_to: str,
    background_tasks: BackgroundTasks,
    db: AsyncSession = Depends(get_db),
    current_admin: UserDB = Depends(get_current_admin_user),
):
    """Send a test email"""
    background_tasks.add_task(send_test_email, email_to=email_to)
    return {"message": f"Test email sent to {email_to} in the background."}
