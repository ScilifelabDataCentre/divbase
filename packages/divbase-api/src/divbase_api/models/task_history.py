"""
Task history DB Model. Summarizes tasks run by Celery without storing all details.
"""

from datetime import datetime
from enum import StrEnum
from typing import TYPE_CHECKING

from sqlalchemy import Column, DateTime, Enum, ForeignKey, Integer, LargeBinary, String, Text
from sqlalchemy.orm import Mapped, mapped_column, relationship

from divbase_api.models.base import Base, BaseDBModel

if TYPE_CHECKING:
    from divbase_api.models.projects import ProjectDB
    from divbase_api.models.users import UserDB


class TaskStatus(StrEnum):
    """
    Helper class that contains the valid Celery task states.
    Used by TaskHistoryDB to set the status column.
    """

    PENDING = "pending"
    STARTED = "started"
    SUCCESS = "success"
    FAILURE = "failure"
    RETRY = "retry"
    REVOKED = "revoked"


class TaskHistoryDB(BaseDBModel):
    """
    DB model for history of tasks executed by Celery.
    """

    __tablename__ = "task_history"

    task_id: Mapped[str] = mapped_column(String, primary_key=True, index=True)
    user_id: Mapped[int] = mapped_column(
        Integer,
        ForeignKey("user.id", ondelete="CASCADE"),
        nullable=True,  # nullable so that cronjob tasks can use user_id None
        index=True,
    )
    project_id: Mapped[int | None] = mapped_column(
        Integer,
        ForeignKey("project.id", ondelete="CASCADE"),
        nullable=True,  # nullable so that cronjob tasks can use project_id None
        index=True,
    )
    status: Mapped[TaskStatus] = mapped_column(Enum(TaskStatus), nullable=False, default=TaskStatus.PENDING)
    error_message: Mapped[str] = mapped_column(String, nullable=True)
    result_message: Mapped[str] = mapped_column(String, nullable=True)

    started_at: Mapped[datetime] = mapped_column(DateTime, nullable=True)
    completed_at: Mapped[datetime] = mapped_column(DateTime, nullable=True)

    user: Mapped["UserDB"] = relationship("UserDB", back_populates="task_history")
    project: Mapped["ProjectDB"] = relationship("ProjectDB", back_populates="task_history")
    celery_meta: Mapped["CeleryTaskMeta"] = relationship(
        "CeleryTaskMeta",
        uselist=False,  # one-to-one relationship: one entry in TaskHistoryDB <-> one entry in CeleryTaskMeta
        viewonly=True,  # Read-only since we don't manage CeleryTaskMeta directly, it is initiated and updated by Celery
    )

    @property
    def runtime_seconds(self) -> float | None:
        """
        Calculate runtime in seconds from started_at and completed_at.
        This is for the admin panel. The deserializer calculates this for the task_history CLI separately
        since property cannot be directly used in the CRUD query.
        """
        if self.started_at and self.completed_at:
            return (self.completed_at - self.started_at).total_seconds()
        return None


class CeleryTaskMeta(Base):
    """
    DB model for the celery_taskmeta table (auto-created by celery in tasks.py in app.conf.update). Not referenced in
    models.__init__py since celery handles table creation.

    Note: this should not inherit from BaseDBModel, since the db table this is referring to is created by celery
    """

    __tablename__ = "celery_taskmeta"

    id = Column(Integer, primary_key=True, autoincrement=True)
    task_id = Column(String(155), ForeignKey("task_history.task_id"), unique=True, index=True)
    status = Column(String(50))
    result = Column(LargeBinary)
    date_done = Column(DateTime)
    traceback = Column(Text)
    name = Column(String(155))
    args = Column(LargeBinary)
    kwargs = Column(LargeBinary)
    worker = Column(String(155))
    retries = Column(Integer)
    queue = Column(String(155))
