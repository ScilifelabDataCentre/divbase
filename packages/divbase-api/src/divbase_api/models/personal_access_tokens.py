"""
Personal Access Token (PAT) DB Model.
"""

from datetime import datetime
from typing import TYPE_CHECKING

from sqlalchemy import DateTime, ForeignKey, String
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import Mapped, mapped_column, relationship

from divbase_api.models.base import BaseDBModel

if TYPE_CHECKING:
    from divbase_api.models.users import UserDB


class PersonalAccessTokenDB(BaseDBModel):
    """
    DB Model for Personal Access Tokens (PATs).

    PATs allow users to authenticate programmatically (e.g., in HPC job submissions/pipeline environments)
    without storing their password. PAT used as a Bearer token.

    permissions: JSONB describing the scopes/permissions of the PAT.
        Structure of JSONB defined by PATPermissions pydantic model in divbase-lib

    id, created_at and updated_at are inherited from BaseDBModel.
    """

    __tablename__ = "personal_access_token"

    user_id: Mapped[int] = mapped_column(ForeignKey("user.id", ondelete="CASCADE"), index=True)
    name: Mapped[str] = mapped_column(String(100))
    description: Mapped[str | None] = mapped_column(String(500), nullable=True)
    hashed_token: Mapped[str] = mapped_column(String(64), index=True, unique=True)  # SHA-256 hex = 64 chars
    permissions: Mapped[dict] = mapped_column(JSONB, nullable=False)

    expires_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True, default=None)
    last_used_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True, default=None)

    is_deleted: Mapped[bool] = mapped_column(default=False)
    date_deleted: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True, default=None)

    user: Mapped["UserDB"] = relationship("UserDB", back_populates="personal_access_tokens")

    def __repr__(self) -> str:
        return f"<PersonalAccessTokenDB id={self.id}, user_id={self.user_id}, name={self.name!r}, is_deleted={self.is_deleted}>"
