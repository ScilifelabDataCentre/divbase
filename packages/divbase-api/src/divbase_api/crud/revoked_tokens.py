"""
CRUD operations for working with revoked JWTs.
"""

import logging

from sqlalchemy import select
from sqlalchemy.exc import IntegrityError
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.models.revoked_tokens import RevokedTokenDB, TokenRevokeReason
from divbase_api.security import TokenType

logger = logging.getLogger(__name__)


async def token_is_revoked(db: AsyncSession, token_jti: str) -> bool:
    """Check if a token has been revoked."""
    stmt = select(RevokedTokenDB).where(RevokedTokenDB.token_jti == token_jti)
    result = await db.execute(stmt)
    return result.scalar_one_or_none() is not None


async def revoke_token_on_logout(db: AsyncSession, token_jti: str, user_id: int) -> None:
    """Revoke a refresh token with reason that user has logged out."""
    await _revoke_token(
        db=db,
        token_jti=token_jti,
        token_type=TokenType.REFRESH,
        user_id=user_id,
        reason=TokenRevokeReason.LOGOUT,
    )


async def revoke_used_password_reset_token(db: AsyncSession, token_jti: str, user_id: int) -> None:
    """Revoke a password reset token with reason that it has been used."""
    await _revoke_token(
        db=db,
        token_jti=token_jti,
        token_type=TokenType.PASSWORD_RESET,
        user_id=user_id,
        reason=TokenRevokeReason.TOKEN_USED,
    )


async def _revoke_token(
    db: AsyncSession, token_jti: str, token_type: TokenType, user_id: int, reason: TokenRevokeReason
) -> None:
    """Helper fun to add a token to the revoked tokens list."""
    logger.info(f"Revoking token jti: {token_jti} of type: {token_type} for user id: {user_id} due to {reason}")
    revoked_token = RevokedTokenDB(
        token_jti=token_jti,
        token_type=token_type,
        user_id=user_id,
        revoked_reason=reason,
    )

    db.add(revoked_token)
    try:
        await db.commit()
    except IntegrityError as e:
        await db.rollback()
        error_details = str(e.orig).lower()
        if "token_jti" in error_details and "unique constraint" in error_details:
            logger.warning(f"Token with jti: {token_jti} was already revoked, ignoring duplicate")
        else:
            logger.error(f"Unexpected integrity error when attempt to revoke token: {e}")
            raise
