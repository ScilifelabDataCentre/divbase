"""
CRUD operations for projects.
"""

from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.models.projects import ProjectDB
from divbase_api.schemas.projects import ProjectCreate


async def create_project(db: AsyncSession, proj_data: ProjectCreate) -> ProjectDB:
    """Create a new project"""
    project = ProjectDB(**proj_data.model_dump())
    db.add(project)
    await db.commit()
    await db.refresh(project)
    return project


async def get_project_by_id(db: AsyncSession, id: int) -> ProjectDB | None:
    """Get project by ID."""
    return await db.get(entity=ProjectDB, ident=id)


async def get_all_projects(db: AsyncSession, limit: int = 1000) -> list[ProjectDB] | list[None]:
    """Get all projects."""
    stmt = select(ProjectDB).limit(limit)
    result = await db.execute(stmt)
    return list(result.scalars().all())
