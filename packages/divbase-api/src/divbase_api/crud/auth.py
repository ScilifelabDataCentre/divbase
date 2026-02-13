"""
Authentication-related CRUD operations.
"""

import logging
from datetime import datetime, timezone

from fastapi import Response
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.revoked_tokens import token_is_revoked
from divbase_api.crud.users import get_user_by_email, get_user_by_id, get_user_by_id_or_raise
from divbase_api.exceptions import AuthenticationError
from divbase_api.models.users import UserDB
from divbase_api.schemas.users import UserPasswordUpdate
from divbase_api.security import TokenType, get_password_hash, verify_password, verify_token

logger = logging.getLogger(__name__)


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


async def verify_user_from_access_token(db: AsyncSession, token: str) -> UserDB | None:
    """
    Attempt to verify a user from their access token, returning the UserDB if successful.

    Returns None if verification fails.
    """
    token_data = verify_token(token=token, desired_token_type=TokenType.ACCESS)
    if not token_data:
        return None
    user = await get_user_by_id(db=db, id=token_data.user_id)
    if not user:
        logger.warning(f"A valid access token was used for a non-existent user with id: {token_data.user_id}.")
        return None
    return user


async def verify_user_from_refresh_token(db: AsyncSession, token: str) -> UserDB | None:
    """
    Attempt to verify a user from their refresh token.
    Returns the UserDB model if successful and None if verification fails.

    As a refresh token is long lived (compared to access tokens), we add extra validation checks here.
    """
    token_data = verify_token(token=token, desired_token_type=TokenType.REFRESH)
    if not token_data:
        return None

    user = await get_user_by_id(db=db, id=token_data.user_id)
    if not user:
        logger.warning(f"A valid refresh token was used for a non-existent user with id: {token_data.user_id}.")
        return None
    if not user_account_valid(user):
        logger.warning(f"Attempt to use refresh token for invalid user account: {user.email} (id: {user.id})")
        return None
    if user.last_password_change and token_data.issued_at < user.last_password_change:
        logger.info(
            f"Refresh token issued before last password change for user: {user.email} (id: {user.id}). Rejecting token."
        )
        return None

    if await token_is_revoked(db=db, token_jti=token_data.jti):
        logger.warning(
            f"Attempt to use revoked refresh token with jti: {token_data.jti} for user id: {token_data.user_id}"
        )
        return None
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
    response.delete_cookie(TokenType.ACCESS.value)
    response.delete_cookie(TokenType.REFRESH.value)
    return response
