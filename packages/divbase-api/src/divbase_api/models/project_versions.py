from datetime import datetime
from typing import TYPE_CHECKING

from sqlalchemy import DateTime, ForeignKey, String, UniqueConstraint
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import Mapped, mapped_column, relationship

from divbase_api.models.base import BaseDBModel

if TYPE_CHECKING:
    from divbase_api.models.projects import ProjectDB


class ProjectVersionDB(BaseDBModel):
    """
    DB Model for project versions.
    An entry here corresponds to the overall state of all files in the project's bucket
    at a given timestamp. Entries are manually created by the users via API.

    The "files" column stores the mapping a JSON object with structure:
    {
        "file_name1": "version_id1",
        "file_name2": "version_id2"
    }

    id, created_at and updated_at are inherited from BaseDBModel.
    """

    __tablename__ = "project_version"

    name: Mapped[str] = mapped_column(String(100), index=True)
    description: Mapped[str | None] = mapped_column(String(500), nullable=True)

    files: Mapped[dict[str, str]] = mapped_column(JSONB, nullable=False)
    is_deleted: Mapped[bool] = mapped_column(default=False)
    date_deleted: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True, default=None)

    # As versioning on project level, can safely delete entry if project deleted
    project_id: Mapped[int] = mapped_column(ForeignKey("project.id", ondelete="CASCADE"), index=True)
    # But we should allow the user who made this to be deleted without deleting the version entry
    user_id: Mapped[int | None] = mapped_column(ForeignKey("user.id", ondelete="SET NULL"))

    project: Mapped["ProjectDB"] = relationship("ProjectDB", back_populates="project_versions")

    __table_args__ = (UniqueConstraint("name", "project_id", name="unique_name_project"),)

    def __repr__(self) -> str:
        return f"<ProjectVersionDB id={self.id}, project_id={self.project_id}, name={self.name}>"
