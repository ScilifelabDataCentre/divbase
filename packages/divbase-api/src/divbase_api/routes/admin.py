"""
Admin only routes for basic DivBase management tasks.

This maybe fully replaced by something like: Starlette-Admin

NOTE: All routes in here should depend on get_current_admin_user.
"""

import logging

from fastapi import APIRouter, Depends, Query, status
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.projects import add_project_member, create_project
from divbase_api.crud.users import (
    create_user,
    deactivate_user,
    get_all_users,
    get_user_by_id_or_raise,
    reactivate_user,
    revert_soft_delete_user,
    soft_delete_user,
)
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
    db: AsyncSession = Depends(get_db),
    current_admin: UserDB = Depends(get_current_admin_user),
):
    """Create a new regular or admin user."""
    new_user = await create_user(db=db, user_data=user_data, is_admin=is_admin)
    logger.info(f"Admin user: {current_admin.email} created a new user: {new_user.email}")
    return new_user


@admin_router.get("/{user_id}", response_model=UserResponse, status_code=status.HTTP_200_OK)
async def get_user_by_id_endpoint(
    user_id: int, db: AsyncSession = Depends(get_db), current_admin: UserDB = Depends(get_current_admin_user)
):
    user = await get_user_by_id_or_raise(db=db, id=user_id)
    return UserResponse.model_validate(user)


@admin_router.get("/", response_model=list[UserResponse], status_code=status.HTTP_200_OK)
async def get_all_users_endpoint(
    limit: int = 1000, db: AsyncSession = Depends(get_db), current_admin: UserDB = Depends(get_current_admin_user)
):
    users = await get_all_users(db, limit=limit)
    return [UserResponse.model_validate(user) for user in users]


@admin_router.patch("/users/{user_id}/deactivate", response_model=UserResponse, status_code=status.HTTP_200_OK)
async def deactivate_user_endpoint(
    user_id: int,
    db: AsyncSession = Depends(get_db),
    current_admin: UserDB = Depends(get_current_admin_user),
):
    """Deactivate a user"""
    user = await deactivate_user(db=db, user_id=user_id)
    logger.info(f"Admin user: {current_admin.email} deactivated user: {user.email}")
    return user


@admin_router.patch("/users/{user_id}/activate", response_model=UserResponse, status_code=status.HTTP_200_OK)
async def activate_user_endpoint(
    user_id: int,
    db: AsyncSession = Depends(get_db),
    current_admin: UserDB = Depends(get_current_admin_user),
):
    """Activate a user"""
    user = await reactivate_user(db=db, user_id=user_id)
    logger.info(f"Admin user: {current_admin.email} activated user: {user.email}")
    return user


@admin_router.patch("/users/{user_id}/soft_delete", response_model=UserResponse, status_code=status.HTTP_200_OK)
async def soft_delete_user_endpoint(
    user_id: int,
    db: AsyncSession = Depends(get_db),
    current_admin: UserDB = Depends(get_current_admin_user),
):
    """Soft delete a user"""
    user = await soft_delete_user(db=db, user_id=user_id)
    logger.info(f"Admin user: {current_admin.email} soft deleted user: {user.email}")
    return user


@admin_router.patch("/users/{user_id}/revert_soft_delete", response_model=UserResponse, status_code=status.HTTP_200_OK)
async def revert_soft_delete_user_endpoint(
    user_id: int,
    db: AsyncSession = Depends(get_db),
    current_admin: UserDB = Depends(get_current_admin_user),
):
    """Revert soft delete of a user"""
    user = await revert_soft_delete_user(db=db, user_id=user_id)
    logger.info(f"Admin user: {current_admin.email} reverted soft delete for user: {user.email}")
    return user


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
    db: AsyncSession = Depends(get_db),
    current_admin: UserDB = Depends(get_current_admin_user),
):
    """Send a test email"""
    send_test_email(email_to=email_to)
    return {"message": f"Test email sent to {email_to}"}
