"""
CRUD operations for users.
"""

from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.models.users import UserDB
from divbase_api.schemas.users import UserCreate
from divbase_api.security import get_password_hash


async def create_user(db: AsyncSession, user_data: UserCreate, is_admin: bool = False) -> UserDB:
    """Create a new user"""
    user_dict = user_data.model_dump(exclude={"password"})
    hashed_password = get_password_hash(user_data.password)

    user = UserDB(**user_dict, hashed_password=hashed_password, is_admin=is_admin)
    db.add(user)
    await db.commit()
    await db.refresh(user)
    return user


async def get_user_by_id(db: AsyncSession, id: int) -> UserDB | None:
    """Get user by ID."""
    return await db.get(entity=UserDB, ident=id)


async def get_user_by_email(db: AsyncSession, email: str) -> UserDB | None:
    """Get user by email."""
    stmt = select(UserDB).where(UserDB.email == email)
    result = await db.execute(stmt)
    return result.scalar_one_or_none()


async def get_all_users(db: AsyncSession, limit: int = 1000, admins_only: bool = False) -> list[UserDB] | list[None]:
    """Get all users."""
    if admins_only:
        stmt = select(UserDB).where(UserDB.is_admin).limit(limit)
    else:
        stmt = select(UserDB).limit(limit)
    result = await db.execute(stmt)
    return list(result.scalars().all())
