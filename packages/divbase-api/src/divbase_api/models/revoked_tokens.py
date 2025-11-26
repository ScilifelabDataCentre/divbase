"""
Revoked Tokens DB Model

Here we store:
- Revoked refresh tokens (e.g. on user logout, password reset)
- Used password reset tokens (to prevent reuse)
"""

from datetime import datetime
from enum import StrEnum

from sqlalchemy import DateTime, ForeignKey, String, func
from sqlalchemy.orm import Mapped, mapped_column, validates

from divbase_api.models.base import BaseDBModel
from divbase_api.security import TokenType


class TokenRevokeReason(StrEnum):
    """
    All reasons why a token may need to be revoked
    """

    LOGOUT = "logout"  # Refresh token reasons
    TOKEN_USED = "token_used"  # Password reset token reasons
    MANUAL_REVOKE = "manual_revoke"  # Admin manually revoked token


class RevokedTokensDB(BaseDBModel):
    """
    DB Model for revoked tokens (can be of type refresh or password reset).

    id, created_at and updated_at are inherited from BaseDBModel.
    """

    __tablename__ = "revoked_tokens"

    token_jti: Mapped[str] = mapped_column(String(36), index=True)  # UUIDv4 string
    token_type: Mapped[TokenType] = mapped_column(String(50), index=True)
    revoked_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), server_default=func.now())
    revoked_reason: Mapped[TokenRevokeReason] = mapped_column(String(100))

    user_id: Mapped[int] = mapped_column(ForeignKey("user.id", ondelete="CASCADE"), index=True)

    @validates("token_type")
    def validate_token_type(self, key, value):
        """Only allow refresh and password reset tokens to be revoked."""
        if value not in {TokenType.REFRESH, TokenType.PASSWORD_RESET}:
            raise ValueError(f"Only refresh and password reset tokens can be revoked, got {value}")
        return value
