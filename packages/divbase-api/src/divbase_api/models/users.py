"""
User DB Model.
"""

from datetime import datetime
from typing import TYPE_CHECKING

from sqlalchemy import DateTime, String, func
from sqlalchemy.orm import Mapped, mapped_column, relationship

from divbase_api.models.base import BaseDBModel

if TYPE_CHECKING:
    from divbase_api.models.projects import ProjectMembershipDB
    from divbase_api.models.revoked_tokens import RevokedTokenDB
    from divbase_api.models.task_history import TaskHistoryDB


class UserDB(BaseDBModel):
    """
    DB Model for a user.

    id, created_at and updated_at are inherited from BaseDBModel.
    """

    __tablename__ = "user"

    name: Mapped[str] = mapped_column(String(100), index=True)
    email: Mapped[str] = mapped_column(String(50), index=True, unique=True)
    hashed_password: Mapped[str] = mapped_column(String(128))  # supports up to SHA-512

    organisation: Mapped[str] = mapped_column(String(200))
    organisation_role: Mapped[str] = mapped_column(String(100))

    is_admin: Mapped[bool] = mapped_column(default=False)
    is_active: Mapped[bool] = mapped_column(default=True)
    is_deleted: Mapped[bool] = mapped_column(default=False)
    email_verified: Mapped[bool] = mapped_column(default=False)

    date_deleted: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True, default=None)
    last_password_change: Mapped[datetime] = mapped_column(DateTime(timezone=True), server_default=func.now())

    project_memberships: Mapped[list["ProjectMembershipDB"]] = relationship(
        "ProjectMembershipDB", back_populates="user"
    )
    task_history: Mapped[list["TaskHistoryDB"]] = relationship("TaskHistoryDB", back_populates="user")
    revoked_tokens: Mapped[list["RevokedTokenDB"]] = relationship("RevokedTokenDB", back_populates="user")

    def __repr__(self) -> str:
        return f"<UserDB id={self.id}, name={self.name}, email={self.email}, is_admin={self.is_admin}>"
