"""
Authentication-related CRUD operations.
"""

from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.api_config import settings
from divbase_api.crud.users import get_user_by_email, get_user_by_id_or_raise
from divbase_api.exceptions import AuthenticationError
from divbase_api.models.users import UserDB
from divbase_api.security import create_email_verification_token, verify_password
from divbase_api.services.email_sender import send_verification_email


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


def send_user_verification_email(user: UserDB) -> None:
    """
    Send a verification email to a user who has just registered to verify their email address.

    TODO: Consider if this should be a background task or sent to celery workers?
    """
    verification_token, _ = create_email_verification_token(subject=user.id)
    verification_url = f"{settings.api.frontend_base_url}/auth/verify-email?token={verification_token}"

    send_verification_email(email_to=user.email, verification_url=verification_url)
