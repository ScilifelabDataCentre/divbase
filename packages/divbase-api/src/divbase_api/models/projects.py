"""
Project DB Model.
"""

from enum import StrEnum

from sqlalchemy import BigInteger
from sqlalchemy.orm import Mapped, mapped_column

from divbase_api.models.base import BaseDBModel


class ProjectRoles(StrEnum):
    """
    Roles assignable to a user on a project.

    - READ: Can view project and download files
    - EDIT: Can add/edit/(soft)delete files + submit queries.
    - MANAGE: Can EDIT + add/remove users and change their roles.

    Note: Admin users do not need to be assigned to projects,
    but have access to all projects.
    """

    READ = "read"
    EDIT = "edit"
    MANAGE = "manage"


class ProjectDB(BaseDBModel):
    """
    DB Model for a project.

    id, created_at and updated_at are inherited from BaseDBModel.

    boto3 lib gives bucket storage used back as integer in bytes,
    hence current choice below to match.
    """

    __tablename__ = "project"

    name: Mapped[str] = mapped_column(index=True, unique=True)
    description: Mapped[str] = mapped_column(nullable=True)
    bucket_name: Mapped[str] = mapped_column(index=True, unique=True)
    is_active: Mapped[bool] = mapped_column(default=True)
    # if not BigInteger, max size would be ~2.1 GB
    storage_quota_bytes: Mapped[int] = mapped_column(BigInteger)
    storage_used_bytes: Mapped[int] = mapped_column(BigInteger, default=0)
