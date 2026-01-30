"""
Task history DB Model. Summarizes tasks run by Celery without storing all details.
"""

from datetime import datetime
from typing import TYPE_CHECKING

from sqlalchemy import Column, DateTime, ForeignKey, Integer, LargeBinary, String, Text
from sqlalchemy.orm import Mapped, mapped_column, relationship

from divbase_api.models.base import Base, BaseDBModel

if TYPE_CHECKING:
    from divbase_api.models.projects import ProjectDB
    from divbase_api.models.users import UserDB


class TaskHistoryDB(BaseDBModel):
    """
    DB model for history of tasks executed by Celery.

    The ID field inherited from BaseDBModel serves as the primary key for this table and is the DivBase job ID.
    The task_id field corresponds to the Celery task ID (UUID string) and is nullable since there are cases such as the bcftools pipe task
    that need the DivBase job ID as input, meaning that a table entry is needed before the Celery task is created.
    """

    __tablename__ = "task_history"

    task_id: Mapped[str | None] = mapped_column(String, index=True, unique=True, nullable=True)
    user_id: Mapped[int | None] = mapped_column(
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

    user: Mapped["UserDB"] = relationship("UserDB", back_populates="task_history")
    project: Mapped["ProjectDB"] = relationship("ProjectDB", back_populates="task_history")
    celery_meta: Mapped["CeleryTaskMeta"] = relationship(
        "CeleryTaskMeta",
        uselist=False,  # one-to-one relationship: one entry in TaskHistoryDB <-> one entry in CeleryTaskMeta
        viewonly=True,  # Read-only since we don't manage CeleryTaskMeta directly, it is initiated and updated by Celery
    )
    started_at_table: Mapped["TaskStartedAtDB"] = relationship(
        "TaskStartedAtDB", primaryjoin="TaskHistoryDB.task_id==foreign(TaskStartedAtDB.task_id)"
    )


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


class TaskStartedAtDB(BaseDBModel):
    """
    DB model for storing the start time of tasks. Not provided by CeleryTaskMeta. Avoids potential race conditions by not having to make multiple updates to entries in TaskHistoryDB.
    """

    __tablename__ = "task_started_at"

    task_id: Mapped[str] = mapped_column(String, index=True, unique=True)
    started_at: Mapped[datetime] = mapped_column(DateTime, nullable=True)
