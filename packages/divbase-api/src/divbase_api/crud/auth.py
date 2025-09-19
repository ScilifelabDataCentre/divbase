"""
Authentication-related CRUD operations.
"""

from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.users import get_user_by_email
from divbase_api.models.users import UserDB
from divbase_api.security import verify_password


async def authenticate_user(db: AsyncSession, email: str, password: str) -> UserDB | None:
    """Authenticate user by email and password when logging in."""
    user = await get_user_by_email(db, email)
    if not user:
        return None

    if not verify_password(plain_password=password, hashed_password=user.hashed_password):
        return None

    return user
