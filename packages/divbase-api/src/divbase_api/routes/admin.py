"""
Admin only routes.

NOTE: All routes in here should depend on get_current_admin_user.
"""

import logging

from fastapi import APIRouter, Depends, HTTPException, Query, status
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.users import create_user, get_user_by_id
from divbase_api.db import get_db
from divbase_api.deps import get_current_admin_user
from divbase_api.models.users import UserDB
from divbase_api.schemas.users import UserCreate, UserResponse

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


@admin_router.patch("/users/{user_id}/deactivate", response_model=UserResponse, status_code=status.HTTP_200_OK)
async def deactivate_user_endpoint(
    user_id: int,
    db: AsyncSession = Depends(get_db),
    current_admin: UserDB = Depends(get_current_admin_user),
):
    """Deactivate a user"""
    user = await get_user_by_id(db=db, id=user_id)
    if not user:
        raise HTTPException(status_code=404, detail="User not found")
    user.is_active = False
    await db.commit()
    await db.refresh(user)
    logger.info(f"Admin user: {current_admin.email} deactivated user: {user.email}")
    return user


@admin_router.patch("/users/{user_id}/activate", response_model=UserResponse, status_code=status.HTTP_200_OK)
async def activate_user_endpoint(
    user_id: int,
    db: AsyncSession = Depends(get_db),
    current_admin: UserDB = Depends(get_current_admin_user),
):
    """Activate a user"""
    user = await get_user_by_id(db=db, id=user_id)
    if not user:
        raise HTTPException(status_code=404, detail="User not found")
    user.is_active = True
    await db.commit()
    await db.refresh(user)
    logger.info(f"Admin user: {current_admin.email} activated user: {user.email}")
    return user
