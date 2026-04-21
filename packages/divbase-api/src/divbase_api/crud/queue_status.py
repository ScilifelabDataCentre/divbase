"""
Get the current status of the queuing system as defined in the QueueStatus model.

This is used to determine if we should reject new tasks from being created when the queue is closed for/ahead of planned maintenance.

CREATE, EDIT and DELETE all handled by admin panel.
"""

from datetime import datetime

from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.exceptions import QueueClosedError
from divbase_api.models.queue_status import QueueStatusDB


async def check_queue_closed_for_new_tasks(db: AsyncSession, is_admin: bool) -> None:
    """
    Checks if the queue is closed for new tasks based on the QueueStatus model.

    Will raise a QueueClosedError if so.
    """
    if is_admin:
        # admins can always submit new tasks to allow for testing and validation during maintenance periods
        return

    result = await db.execute(select(QueueStatusDB).filter_by(id=1))
    status = result.scalar_one_or_none()

    if not status or not status.is_closed:
        return

    if not status.scheduled_start:  # effective immediately if no scheduled start time
        raise QueueClosedError(message=status.reason_for_users)

    # if the queue is scheduled to be closed in the future, we can still allow new tasks to be created until that time
    if status.scheduled_start < datetime.now(tz=status.scheduled_start.tzinfo):
        raise QueueClosedError(message=status.reason_for_users)
