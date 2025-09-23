"""
Frontend routes for user management.

These routes will return Template Responses.
"""

from fastapi import APIRouter, Depends, status
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.users import get_user_by_id
from divbase_api.db import get_db
from divbase_api.schemas.users import UserResponse

frontend_users_router = APIRouter()


@frontend_users_router.get("/{user_id}", response_model=UserResponse, status_code=status.HTTP_200_OK)
async def get_user_by_id_endpoint(user_id: int, db: AsyncSession = Depends(get_db)):
    users = await get_user_by_id(db=db, id=user_id)
    return users
