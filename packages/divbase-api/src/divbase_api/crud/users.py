"""
CRUD operations for users.
"""

from datetime import datetime, timezone

from fastapi import HTTPException
from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.exceptions import UserRegistrationError
from divbase_api.models.users import UserDB
from divbase_api.schemas.users import UserCreate, UserUpdate
from divbase_api.security import get_password_hash


async def get_user_by_id(db: AsyncSession, id: int) -> UserDB | None:
    """Get user by ID."""
    return await db.get(entity=UserDB, ident=id)


async def get_user_by_id_or_raise(db: AsyncSession, id: int) -> UserDB:
    """Get user by ID, but raise 404 if user not found."""
    user = await db.get(entity=UserDB, ident=id)
    if not user:
        raise HTTPException(status_code=404, detail="User not found")
    return user


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


async def create_user(
    db: AsyncSession, user_data: UserCreate, is_admin: bool = False, email_verified: bool = False
) -> UserDB:
    """Create a new user"""
    proposed_email = user_data.email.lower()
    current_user = await get_user_by_email(db=db, email=proposed_email)
    if current_user:
        # TODO this could be changed to say, you'll recieve an email to verify account
        # but we actually send a your email is already verified email
        # To prevent email enumeration attacks
        raise UserRegistrationError(
            internal_logging_message=f"Attempt made to register new account with existing email: {proposed_email}"
        )

    user_dict = user_data.model_dump(exclude={"password"})
    hashed_password = get_password_hash(user_data.password)

    user = UserDB(**user_dict, hashed_password=hashed_password, is_admin=is_admin, email_verified=email_verified)
    db.add(user)
    await db.commit()
    await db.refresh(user)
    return user


async def deactivate_user(db: AsyncSession, user_id: int) -> UserDB:
    """Deactivate a user"""
    user = await get_user_by_id_or_raise(db=db, id=user_id)
    user.is_active = False
    await db.commit()
    await db.refresh(user)
    return user


async def reactivate_user(db: AsyncSession, user_id: int) -> UserDB:
    """Reactivate a user"""
    user = await get_user_by_id_or_raise(db=db, id=user_id)
    user.is_active = True
    await db.commit()
    await db.refresh(user)
    return user


async def soft_delete_user(db: AsyncSession, user_id: int) -> UserDB:
    """Soft delete a user"""
    user = await get_user_by_id_or_raise(db=db, id=user_id)
    user.is_deleted = True
    user.is_active = False  # also deactivate as deleted users should not be active
    user.date_deleted = datetime.now(tz=timezone.utc)
    await db.commit()
    await db.refresh(user)
    return user


async def revert_soft_delete_user(db: AsyncSession, user_id: int) -> UserDB:
    """Revert soft delete of a user"""
    user = await get_user_by_id_or_raise(db=db, id=user_id)
    user.is_deleted = False
    user.is_active = True  # reactivate the user if deactivated too.
    user.date_deleted = None
    await db.commit()
    await db.refresh(user)
    return user


async def update_user_profile(db: AsyncSession, user_data: UserUpdate, user_id: int) -> UserDB:
    """Used by a regular user to update their own profile information."""
    user = await get_user_by_id_or_raise(db=db, id=user_id)

    update_data = user_data.model_dump(exclude_unset=True)
    for field, value in update_data.items():
        setattr(user, field, value)

    await db.commit()
    await db.refresh(user)
    return user
