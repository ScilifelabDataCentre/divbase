"""
Routes for user management.
"""

from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.users import get_all_users, get_user_by_id
from divbase_api.db import get_db
from divbase_api.schemas.users import UserResponse

users_router = APIRouter()


@users_router.get("/{user_id}", response_model=UserResponse, status_code=status.HTTP_200_OK)
async def get_user_by_id_endpoint(user_id: int, db: AsyncSession = Depends(get_db)):
    user = await get_user_by_id(db=db, id=user_id)
    if not user:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="User not found")
    return UserResponse.model_validate(user)


@users_router.get("/", response_model=list[UserResponse], status_code=status.HTTP_200_OK)
async def get_all_users_endpoint(db: AsyncSession = Depends(get_db), limit: int = 1000):
    users = await get_all_users(db, limit=limit)
    return [UserResponse.model_validate(user) for user in users]
