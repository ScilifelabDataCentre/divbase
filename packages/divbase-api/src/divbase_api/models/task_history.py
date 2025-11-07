"""
Task history DB Model. Summarizes tasks run by Celery without storing all details.
"""

import enum
from typing import TYPE_CHECKING

from sqlalchemy import Enum, ForeignKey, Integer, String
from sqlalchemy.orm import Mapped, mapped_column, relationship

from divbase_api.models.base import BaseDBModel

if TYPE_CHECKING:
    from divbase_api.models.projects import ProjectDB
    from divbase_api.models.users import UserDB


class TaskStatus(enum.Enum):
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
    user_id: Mapped[int] = mapped_column(Integer, ForeignKey("user.id", ondelete="CASCADE"), nullable=False, index=True)
    project_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("project.id", ondelete="CASCADE"), nullable=False, index=True
    )
    status: Mapped[TaskStatus] = mapped_column(Enum(TaskStatus), nullable=False, default=TaskStatus.PENDING)
    error_message: Mapped[str] = mapped_column(String, nullable=True)
    result_message: Mapped[str] = mapped_column(String, nullable=True)

    user: Mapped["UserDB"] = relationship("UserDB", back_populates="task_history")
    project: Mapped["ProjectDB"] = relationship("ProjectDB", back_populates="task_history")
