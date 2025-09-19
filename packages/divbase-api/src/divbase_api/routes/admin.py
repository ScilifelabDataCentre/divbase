"""
Admin only routes.

TODO: These routes will need a "get_current_admin_user" dependency.

"""

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.users import create_user, get_user_by_id
from divbase_api.db import get_db
from divbase_api.schemas.users import UserCreate, UserResponse

admin_router = APIRouter()


@admin_router.post("/users/", response_model=UserResponse)
async def create_user_endpoint(
    user_data: UserCreate,
    is_admin: bool = Query(False, description="Set to true to create an admin user"),
    db: AsyncSession = Depends(get_db),
):
    """Create a new regular or admin user."""
    new_user = await create_user(db=db, user_data=user_data, is_admin=is_admin)
    return new_user


@admin_router.patch("/users/{user_id}/deactivate", response_model=UserResponse)
async def deactivate_user_endpoint(
    user_id: int,
    db: AsyncSession = Depends(get_db),
):
    """Deactivate a user"""
    user = await get_user_by_id(db=db, id=user_id)
    if not user:
        raise HTTPException(status_code=404, detail="User not found")
    user.is_active = False
    await db.commit()
    await db.refresh(user)
    return user
