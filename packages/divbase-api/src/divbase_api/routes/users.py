"""
Routes for user management.
"""

from fastapi import APIRouter, Depends, status
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.users import create_user, get_all_users, get_user_by_id
from divbase_api.db import get_db
from divbase_api.schemas.users import UserCreate, UserResponse

users_router = APIRouter()


@users_router.get("/{user_id}", response_model=UserResponse, status_code=status.HTTP_200_OK)
async def get_user_by_id_endpoint(user_id: int, db: AsyncSession = Depends(get_db)):
    users = await get_user_by_id(db=db, id=user_id)
    return users


@users_router.get("/", response_model=list[UserResponse], status_code=status.HTTP_200_OK)
async def get_all_users_endpoint(db: AsyncSession = Depends(get_db), limit: int = 1000):
    users = await get_all_users(db, limit=limit)
    return users


@users_router.post("/", response_model=UserResponse, status_code=status.HTTP_201_CREATED)
async def create_user_endpoint(user_data: UserCreate, db: AsyncSession = Depends(get_db)):
    new_user = await create_user(db=db, user_data=user_data)
    return new_user
