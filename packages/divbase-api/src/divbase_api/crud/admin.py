"""
CRUD operations for admin functionalities.
"""

from sqlalchemy import func, select
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.models.projects import ProjectDB
from divbase_api.models.users import UserDB
from divbase_api.schemas.admin import SystemStatsResponse


async def get_system_stats(db: AsyncSession) -> SystemStatsResponse:
    """Get system statistics for admin dashboard."""

    total_users_stmt = select(func.count(UserDB.id))
    total_users = await db.execute(total_users_stmt)

    inactive_users_stmt = select(func.count(UserDB.id)).filter(~UserDB.is_active)
    inactive_users = await db.execute(inactive_users_stmt)

    admin_users_stmt = select(func.count(UserDB.id)).where(UserDB.is_admin)
    admin_users = await db.execute(admin_users_stmt)

    total_projects_stmt = select(func.count(ProjectDB.id))
    total_projects = await db.execute(total_projects_stmt)

    active_projects_stmt = select(func.count(ProjectDB.id)).filter(ProjectDB.is_active)
    active_projects = await db.execute(active_projects_stmt)

    return SystemStatsResponse(
        total_users=total_users.scalar() or 0,
        inactive_users=inactive_users.scalar() or 0,
        admin_users=admin_users.scalar() or 0,
        total_projects=total_projects.scalar() or 0,
        active_projects=active_projects.scalar() or 0,
        divbase_version="TODO",
        last_updated="TODO",
        environment="TODO",
    )
