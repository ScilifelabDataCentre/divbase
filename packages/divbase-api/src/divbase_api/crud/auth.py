"""
Authentication-related CRUD operations.
"""

from datetime import datetime, timezone

from fastapi import Response
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.users import get_user_by_email, get_user_by_id_or_raise
from divbase_api.exceptions import AuthenticationError
from divbase_api.models.users import UserDB
from divbase_api.schemas.users import UserPasswordUpdate
from divbase_api.security import TokenType, get_password_hash, verify_password


async def authenticate_user(db: AsyncSession, email: str, password: str) -> UserDB:
    """
    Authenticate user by email and password when logging in and validates user
    has access to the system (not deleted, active, email verified).

    Raises AuthenticationError if authentication fails.
    """
    generic_error_msg = "Invalid email or password or user account does not exist."
    user = await get_user_by_email(db, email)

    if not user:
        raise AuthenticationError(message=generic_error_msg)

    if not verify_password(plain_password=password, hashed_password=user.hashed_password):
        raise AuthenticationError(message=generic_error_msg)

    if user.is_deleted or not user.is_active:
        raise AuthenticationError(message=generic_error_msg)

    if not user.email_verified:
        raise AuthenticationError(
            message="Email address not verified, check your inbox or visit the DivBase website to resend a new verification email."
        )

    return user


def user_account_valid(user: UserDB) -> bool:
    """Check if user account is valid (active, not deleted, email verified)."""
    return user.is_active and not user.is_deleted and user.email_verified


async def check_user_email_verified(db: AsyncSession, id: int) -> bool:
    """Check if a user's email is verified."""
    user = await get_user_by_id_or_raise(db=db, id=id)
    return user.email_verified


async def confirm_user_email(db: AsyncSession, id: int) -> UserDB:
    """Update user to set email_verified to True."""
    user = await get_user_by_id_or_raise(db=db, id=id)
    user.email_verified = True
    await db.commit()
    await db.refresh(user)
    return user


async def update_user_password(db: AsyncSession, user_id: int, password_data: UserPasswordUpdate) -> UserDB:
    """change the user's password."""
    user = await get_user_by_id_or_raise(db=db, id=user_id)

    hashed_password = get_password_hash(password_data.password)
    user.hashed_password = hashed_password
    user.last_password_change = datetime.now(tz=timezone.utc)
    await db.commit()
    await db.refresh(user)
    return user


def delete_auth_cookies(response: Response) -> Response:
    """Helper to delete auth cookies from a response (e.g. on logout)."""
    # TODO - when token blacklisting is implemented, blacklist the tokens here too. (make async at that point too...)
    response.delete_cookie(TokenType.ACCESS.value)
    response.delete_cookie(TokenType.REFRESH.value)
    return response
