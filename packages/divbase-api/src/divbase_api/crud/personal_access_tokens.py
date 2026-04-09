"""
CRUD operations for Personal Access Tokens (PATs).
"""

import logging
from datetime import datetime, timezone

from pydantic import SecretStr
from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.auth import user_account_valid
from divbase_api.crud.users import get_user_by_id
from divbase_api.models.personal_access_tokens import PersonalAccessTokenDB
from divbase_api.models.users import UserDB
from divbase_api.security import generate_personal_access_token, hash_personal_access_token

logger = logging.getLogger(__name__)


async def create_personal_access_token(
    db: AsyncSession,
    user_id: int,
    name: str,
    description: str | None = None,
    permissions: dict | None = None,
    expires_at: datetime | None = None,
) -> tuple[PersonalAccessTokenDB, SecretStr]:
    """
    Create a new personal access token for a user.
    Stores a hashed version of the token in the database,
    and returns the plaintext token for the user to copy (only available at creation time).
    """
    raw_pat = generate_personal_access_token()
    hashed_pat = hash_personal_access_token(raw_pat)

    pat = PersonalAccessTokenDB(
        user_id=user_id,
        name=name,
        description=description,
        hashed_token=hashed_pat,
        permissions=permissions,
        expires_at=expires_at,
    )
    db.add(pat)
    await db.commit()
    await db.refresh(pat)
    return pat, raw_pat


async def verify_user_from_personal_access_token(
    db: AsyncSession, raw_pat: SecretStr
) -> tuple[UserDB, PersonalAccessTokenDB] | None:
    """
    Verify a PAT. Includes validation of the user account/db entry

    Returns the owning user and PAT DB entry if valid.
    Returns None if invalid.
    """
    hashed = hash_personal_access_token(raw_pat)

    stmt = select(PersonalAccessTokenDB).where(PersonalAccessTokenDB.hashed_token == hashed)
    result = await db.execute(stmt)
    hashed_pat = result.scalar_one_or_none()

    if not hashed_pat:
        return None

    if not hashed_pat.is_deleted:
        logger.info(f"PAT id={hashed_pat.id} rejected, as it is soft deleted.")
        return None

    if hashed_pat.expires_at and hashed_pat.expires_at < datetime.now(timezone.utc):
        logger.info(f"PAT id={hashed_pat.id} rejected, as has expired.")
        return None

    user = await get_user_by_id(db=db, id=hashed_pat.user_id)
    if not user:
        logger.warning(f"PAT id={hashed_pat.id} rejected, it references a non-existent user id={hashed_pat.user_id}.")
        return None

    if not user_account_valid(user):
        logger.info(f"PAT id={hashed_pat.id} rejected, as user account is invalid: {user.email}")
        return None

    hashed_pat.last_used_at = datetime.now(timezone.utc)
    await db.commit()

    return user, hashed_pat


async def get_users_personal_access_tokens(db: AsyncSession, user_id: int) -> list[PersonalAccessTokenDB]:
    """
    Get all personal access tokens belonging to the given user.

    (soft deleted PATs are not included.)
    """
    stmt = select(PersonalAccessTokenDB).where(
        PersonalAccessTokenDB.user_id == user_id,
        PersonalAccessTokenDB.is_deleted == False,  # noqa: E712
    )
    result = await db.execute(stmt)
    return list(result.scalars().all())


async def soft_delete_personal_access_token(db: AsyncSession, pat_id: int, user_id: int) -> bool:
    """
    Soft delete a Personal Access Token by ID, scoped to a user (fn to be used in user-facing operations).
    Admins can soft and hard delete any PAT via the admin_panel instead.

    Returns False if the PAT does not exist or does not belong to the user.
    True if successfully soft deleted.
    """
    stmt = select(PersonalAccessTokenDB).where(
        PersonalAccessTokenDB.id == pat_id,
        PersonalAccessTokenDB.user_id == user_id,
    )
    result = await db.execute(stmt)
    pat = result.scalar_one_or_none()
    if not pat:
        return False

    pat.is_deleted = True
    pat.date_deleted = datetime.now(tz=timezone.utc)
    await db.commit()
    await db.refresh(pat)
    logger.info(f"PAT id={pat.id} soft deleted for user_id={user_id}.")
    return True
