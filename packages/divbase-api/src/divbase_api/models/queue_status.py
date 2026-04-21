"""
QueueStatus db model is used to define whether the DivBase queuing system is closed to new tasks.

The use case for this is prior to a scheduled upgrade, we may want to empty the queue:
1. (Day before upgrade): Go to admin panel, set is_closed to True (optionally schedule when the queue will be closed).
2. (Day of upgrade): Validate queue is empty, perform upgrade,
3. Once happy, go to admin panel, set is_closed to False to allow new tasks to be created again.

NOTE: Admins can submit jobs even when the queue is closed to allow for validation and testing during the maintenance period
"""

from datetime import datetime

from sqlalchemy import CheckConstraint, DateTime, String
from sqlalchemy.orm import Mapped, mapped_column

from divbase_api.models.base import BaseDBModel


class QueueStatusDB(BaseDBModel):
    __tablename__ = "queue_status"

    is_closed: Mapped[bool] = mapped_column(default=False)
    scheduled_start: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)
    reason_for_users: Mapped[str] = mapped_column(String(500))

    # this arg makes the table a singleton, so only 1 row (with id = 1) in the table
    __table_args__ = (CheckConstraint("id = 1", name="queue_status_singleton_row"),)
