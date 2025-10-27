"""
Authentication-related CRUD operations.
"""

from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.config import settings
from divbase_api.crud.users import get_user_by_email, get_user_by_id_or_raise
from divbase_api.models.users import UserDB
from divbase_api.security import create_email_verification_token, verify_password
from divbase_api.services.email_sender import send_verification_email


async def authenticate_user(db: AsyncSession, email: str, password: str) -> UserDB | None:
    """Authenticate user by email and password when logging in."""
    user = await get_user_by_email(db, email)
    if not user or not user.is_active:
        return None

    if not verify_password(plain_password=password, hashed_password=user.hashed_password):
        return None

    return user


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
