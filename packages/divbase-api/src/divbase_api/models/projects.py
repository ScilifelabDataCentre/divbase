"""
Project DB Model.
"""

from enum import StrEnum
from typing import TYPE_CHECKING

from sqlalchemy import BigInteger, ForeignKey, String, UniqueConstraint
from sqlalchemy.orm import Mapped, mapped_column, relationship

from divbase_api.models.base import BaseDBModel

if TYPE_CHECKING:
    from divbase_api.models.users import UserDB


class ProjectDB(BaseDBModel):
    """
    DB Model for a project.

    id, created_at and updated_at are inherited from BaseDBModel.

    boto3 lib gives bucket storage used back as integer in bytes,
    hence current choice below to match.
    """

    __tablename__ = "project"

    name: Mapped[str] = mapped_column(String(100), index=True, unique=True)
    description: Mapped[str] = mapped_column(String(1000), nullable=True)
    bucket_name: Mapped[str] = mapped_column(String(63), index=True, unique=True)
    is_active: Mapped[bool] = mapped_column(default=True)
    # if not BigInteger, max size would be ~2.1 GB
    storage_quota_bytes: Mapped[int] = mapped_column(BigInteger)
    storage_used_bytes: Mapped[int] = mapped_column(BigInteger, default=0)

    memberships: Mapped[list["ProjectMembershipDB"]] = relationship("ProjectMembershipDB", back_populates="project")

    def __repr__(self) -> str:
        return f"<ProjectDB id={self.id}, name={self.name}, bucket_name={self.bucket_name}, is_active={self.is_active}>"


class ProjectRoles(StrEnum):
    """
    Roles assignable to a user in a project.

    - READ: Can view project and download files
    - EDIT: Can add/edit/(soft)delete files + submit queries.
    - MANAGE: Can EDIT + add/remove users and change their roles.

    Note: Admin users do not need to be assigned to projects,
    but have access to all projects.
    """

    READ = "read"
    EDIT = "edit"
    MANAGE = "manage"


class ProjectMembershipDB(BaseDBModel):
    """
    DB Model for the many-to-many relationship between users and projects.

    id, created_at and updated_at are inherited from BaseDBModel.
    """

    __tablename__ = "project_membership"

    # ondelete=cascade ensures if user or project is deleted, their memberships are also deleted.
    user_id: Mapped[int] = mapped_column(ForeignKey("user.id", ondelete="CASCADE"), index=True)
    project_id: Mapped[int] = mapped_column(ForeignKey("project.id", ondelete="CASCADE"), index=True)

    role: Mapped[ProjectRoles] = mapped_column(index=True)

    user: Mapped["UserDB"] = relationship("UserDB", back_populates="project_memberships")
    project: Mapped["ProjectDB"] = relationship("ProjectDB", back_populates="memberships")

    __table_args__ = (UniqueConstraint("user_id", "project_id", name="unique_user_project"),)

    def __repr__(self) -> str:
        return f"<ProjectMembershipDB id={self.id}, user_id={self.user_id}, project_id={self.project_id}, role={self.role}>"
