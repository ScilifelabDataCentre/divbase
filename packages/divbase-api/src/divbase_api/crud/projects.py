"""
CRUD operations for projects.
"""

from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.orm import selectinload

from divbase_api.crud.users import get_user_by_id
from divbase_api.models.projects import ProjectDB, ProjectMembershipDB, ProjectRoles
from divbase_api.schemas.projects import ProjectCreate, UserProjectResponse


async def create_project(db: AsyncSession, proj_data: ProjectCreate) -> ProjectDB:
    """Create a new project"""
    # TODO - add validation, e.g. check unique name, bucket not already being used etc...
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


async def get_user_projects_with_roles(db: AsyncSession, user_id: int) -> list[tuple[ProjectDB, ProjectRoles]]:
    """Get all projects that a user is a member of and their role in that project."""
    stmt = (
        select(ProjectDB, ProjectMembershipDB.role)
        .join(ProjectMembershipDB)
        .where(ProjectMembershipDB.user_id == user_id)
        .options(selectinload(ProjectDB.memberships))
    )
    result = await db.execute(stmt)

    return [(row[0], ProjectRoles(row[1])) for row in result.all()]


async def get_project_with_user_role(
    db: AsyncSession, project_id: int, user_id: int
) -> tuple[ProjectDB | None, ProjectRoles | None]:
    """Get project by ID and the user's role in that project."""
    stmt = (
        select(ProjectDB, ProjectMembershipDB.role)
        .join(ProjectMembershipDB)
        .where(ProjectDB.id == project_id)
        .where(ProjectMembershipDB.user_id == user_id)
    )
    result = await db.execute(stmt)
    row = result.first()

    if row:
        return row[0], ProjectRoles(row[1])
    return None, None


async def add_project_member(
    db: AsyncSession, project_id: int, user_id: int, role: ProjectRoles
) -> ProjectMembershipDB:
    """Add a user to a project with a defined ProjectRole."""
    # TODO - custom exception for these valueerrors.

    project = await get_project_by_id(db=db, id=project_id)
    if not project:
        raise ValueError("Project not found")

    user = await get_user_by_id(db=db, id=user_id)
    if not user:
        raise ValueError("User not found")

    # TODO - check if membership already exists
    membership = ProjectMembershipDB(project_id=project_id, user_id=user_id, role=role)

    db.add(membership)
    await db.commit()
    await db.refresh(membership)
    return membership


def create_user_project_responses(
    user_projects: list[tuple[ProjectDB, ProjectRoles]],
) -> list[UserProjectResponse]:
    """Helper to create user project responses from list of (ProjectDB, ProjectRoles) tuples."""
    project_responses = []
    for project, user_role in user_projects:
        project_response = UserProjectResponse(
            name=project.name,
            description=project.description,
            bucket_name=project.bucket_name,
            id=project.id,
            is_active=project.is_active,
            user_role=user_role,
            storage_quota_bytes=project.storage_quota_bytes,
            storage_used_bytes=project.storage_used_bytes,
        )
        project_responses.append(project_response)
    return project_responses
