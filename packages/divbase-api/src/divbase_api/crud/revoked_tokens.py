"""
CRUD operations for working with revoked JWTs.
"""

from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.models.revoked_tokens import RevokedTokensDB, TokenRevokeReason
from divbase_api.security import TokenType


async def token_is_revoked(db: AsyncSession, token_jti: str) -> bool:
    """Check if a token has been revoked."""
    stmt = select(RevokedTokensDB).where(RevokedTokensDB.token_jti == token_jti)
    result = await db.execute(stmt)
    return result.scalar_one_or_none() is not None


async def revoke_token_on_logout(db: AsyncSession, token_jti: str, user_id: int) -> RevokedTokensDB:
    """Revoke a refresh token with reason that user has logged out."""
    return await _revoke_token(
        db=db,
        token_jti=token_jti,
        token_type=TokenType.REFRESH,
        user_id=user_id,
        reason=TokenRevokeReason.LOGOUT,
    )


async def revoke_used_password_reset_token(db: AsyncSession, token_jti: str, user_id: int) -> RevokedTokensDB:
    """Revoke a password reset token with reason that it has been used."""
    return await _revoke_token(
        db=db,
        token_jti=token_jti,
        token_type=TokenType.PASSWORD_RESET,
        user_id=user_id,
        reason=TokenRevokeReason.TOKEN_USED,
    )


async def _revoke_token(
    db: AsyncSession, token_jti: str, token_type: TokenType, user_id: int, reason: TokenRevokeReason
) -> RevokedTokensDB:
    """Helper fun to add a token to the revoked tokens list."""
    revoked_token = RevokedTokensDB(
        token_jti=token_jti,
        token_type=token_type,
        user_id=user_id,
        revoked_reason=reason,
    )
    db.add(revoked_token)
    await db.commit()
    await db.refresh(revoked_token)
    return revoked_token
