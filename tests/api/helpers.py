"""Shared helpers for DivBase API tests."""

from datetime import datetime, timedelta, timezone
from uuid import uuid4

from pydantic import SecretStr
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.users import create_user, get_user_by_id_or_raise
from divbase_api.schemas.users import UserCreate, UserResponse

TEST_PASSWORD = "badpassword"


async def create_test_user(
    db_session: AsyncSession,
    email: str | None = None,
    password: str = TEST_PASSWORD,
    email_verified: bool = True,
    is_admin: bool = False,
    is_active: bool = True,
    is_soft_deleted: bool = False,
) -> UserResponse:
    """
    Creates a user for testing.

    Last password change is set to be one hour ago. This is because
    verify_user_from_refresh_token() rejects any refresh token issued before last password change (seconds precision).
    """
    if email is None:
        email = f"{uuid4()}@divbase.com"

    user_data = UserCreate(
        name="Test User",
        email=email,
        organisation="DivBase University",
        organisation_role="Developer",
        password=SecretStr(password),
        confirm_password=SecretStr(password),
    )
    user = await create_user(db=db_session, user_data=user_data, is_admin=is_admin, email_verified=email_verified)

    user_row = await get_user_by_id_or_raise(db=db_session, id=user.id)
    user_row.last_password_change = datetime.now(timezone.utc) - timedelta(hours=1)
    if not is_active:
        user_row.is_active = False
    if is_soft_deleted:
        user_row.is_deleted = True
        user_row.date_deleted = datetime.now(timezone.utc)

    await db_session.commit()
    return user
