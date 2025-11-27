"""
Revoked Tokens DB Model

Here we store:
- Revoked refresh tokens (e.g. on user logout, password reset)
- Used password reset tokens (to prevent reuse)
"""

from datetime import datetime
from enum import StrEnum
from typing import TYPE_CHECKING

from sqlalchemy import DateTime, ForeignKey, String, func
from sqlalchemy.orm import Mapped, mapped_column, relationship, validates

from divbase_api.models.base import BaseDBModel
from divbase_api.security import TokenType

if TYPE_CHECKING:
    from divbase_api.models.users import UserDB


class TokenRevokeReason(StrEnum):
    """
    All reasons why a token may need to be revoked

    NOTE: password change or account deactivated/soft deleted do not require
    explicit revocation as we check those conditions (on user db model) when verifying tokens.
    """

    LOGOUT = "logout"  # Refresh token reason
    TOKEN_USED = "token_used"  # Password reset token reason
    MANUAL_REVOKE = "manual_revoke"  # Admin manually revoked token


class RevokedTokenDB(BaseDBModel):
    """
    DB Model for revoked tokens (used during refresh token and password reset flows).

    id, created_at and updated_at are inherited from BaseDBModel.
    """

    __tablename__ = "revoked_token"

    token_jti: Mapped[str] = mapped_column(String(36), index=True, unique=True)  # UUIDv4 string
    token_type: Mapped[TokenType] = mapped_column(index=True)
    revoked_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), server_default=func.now())
    revoked_reason: Mapped[TokenRevokeReason] = mapped_column()

    user_id: Mapped[int | None] = mapped_column(ForeignKey("user.id", ondelete="CASCADE"), index=True, nullable=True)
    user: Mapped["UserDB | None"] = relationship("UserDB", back_populates="revoked_tokens")

    def __repr__(self) -> str:
        return f"<RevokedTokenDB id={self.id}, token_jti={self.token_jti}, token_type={self.token_type}, revoked_reason={self.revoked_reason}>"

    @validates("token_type")
    def validate_token_type(self, key, value):
        """
        Only allow refresh and password reset tokens to be revoked.

        NOTE: Even if an access or email_verify token is added here, it would not have any effect
        as the is no check against this db when those are used.
        """
        if value not in {TokenType.REFRESH, TokenType.PASSWORD_RESET}:
            raise ValueError(f"Only refresh and password reset tokens can be revoked, got {value}")
        return value
