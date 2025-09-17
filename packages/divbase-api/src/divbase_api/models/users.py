"""
User DB Model.
"""

from typing import TYPE_CHECKING

from sqlalchemy import String
from sqlalchemy.orm import Mapped, mapped_column, relationship

from divbase_api.models.base import BaseDBModel

if TYPE_CHECKING:
    from divbase_api.models.projects import ProjectMembershipDB


class UserDB(BaseDBModel):
    """
    DB Model for a user.

    id, created_at and updated_at are inherited from BaseDBModel.
    """

    __tablename__ = "user"

    name: Mapped[str] = mapped_column(String(100), index=True, unique=True)
    email: Mapped[str] = mapped_column(String(50), index=True, unique=True)
    hashed_password: Mapped[str] = mapped_column(String(128))  # supports up to SHA-512

    is_admin: Mapped[bool] = mapped_column(default=False)
    is_active: Mapped[bool] = mapped_column(default=True)

    project_memberships: Mapped[list["ProjectMembershipDB"]] = relationship(
        "ProjectMembershipDB", back_populates="user"
    )

    def __repr__(self) -> str:
        return f"<UserDB id={self.id}, name={self.name}, email={self.email}, is_admin={self.is_admin}>"
