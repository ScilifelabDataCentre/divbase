"""
CRUD operations for users.
"""

from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.models.users import UserDB
from divbase_api.schemas.users import UserCreate


async def create_user(db: AsyncSession, user_data: UserCreate) -> UserDB:
    """Create a new user"""
    user_dict = user_data.model_dump(exclude={"password"})
    hashed_password = user_data.password.get_secret_value() + "_notreallyhashed"  # TODO: Replace with real hashing

    user = UserDB(**user_dict, hashed_password=hashed_password, is_admin=False)
    db.add(user)
    await db.commit()
    await db.refresh(user)
    return user


async def get_user_by_id(db: AsyncSession, id: int) -> UserDB | None:
    """Get user by ID."""
    return await db.get(entity=UserDB, ident=id)


async def get_all_users(db: AsyncSession, limit: int = 1000) -> list[UserDB] | list[None]:
    """Get all users."""
    stmt = select(UserDB).limit(limit)
    result = await db.execute(stmt)
    return list(result.scalars().all())
