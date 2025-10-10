"""
CRUD operations for projects.
"""

from sqlalchemy import exists, select
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.orm import selectinload

from divbase_api.crud.users import get_user_by_id
from divbase_api.exceptions import ProjectCreationError, ProjectMemberNotFoundError, ProjectNotFoundError
from divbase_api.models.projects import ProjectDB, ProjectMembershipDB, ProjectRoles
from divbase_api.schemas.projects import ProjectCreate, ProjectMemberResponse, UserProjectResponse


async def create_project(db: AsyncSession, proj_data: ProjectCreate) -> ProjectDB:
    """Create a new project"""
    await _validate_project_create(db=db, name=proj_data.name, bucket_name=proj_data.bucket_name)

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


async def get_project_with_user_role(db: AsyncSession, project_id: int, user_id: int) -> tuple[ProjectDB, ProjectRoles]:
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
    raise ProjectNotFoundError("Project not found or the user has no access")


async def get_project_id_from_name(db: AsyncSession, project_name: str) -> int | None:
    """Get a project's ID from its name."""
    stmt = select(ProjectDB.id).where(ProjectDB.name == project_name)
    result = await db.execute(stmt)
    project_id = result.scalar_one_or_none()
    return project_id


async def add_project_member(
    db: AsyncSession, project_id: int, user_id: int, role: ProjectRoles
) -> ProjectMembershipDB:
    """Add a user to a project with a defined ProjectRole."""
    # TODO - custom exception for these valueerrors.

    project = await get_project_by_id(db=db, id=project_id)
    if not project:
        raise ProjectNotFoundError()

    user = await get_user_by_id(db=db, id=user_id)
    if not user:
        raise ValueError("User not found")

    already_member = await is_user_member_of_project(db=db, project_id=project_id, user_id=user_id)
    if already_member:
        raise ValueError("User is already a member of this project")

    membership = ProjectMembershipDB(project_id=project_id, user_id=user_id, role=role)

    db.add(membership)
    await db.commit()
    await db.refresh(membership)
    return membership


async def is_user_member_of_project(db: AsyncSession, project_id: int, user_id: int) -> bool:
    """Check if a user is already a member of a project."""
    stmt = select(ProjectMembershipDB).where(
        ProjectMembershipDB.project_id == project_id,
        ProjectMembershipDB.user_id == user_id,
    )
    result = await db.execute(stmt)
    membership = result.scalar_one_or_none()
    return membership is not None


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


async def get_project_members(db: AsyncSession, project_id: int) -> list[ProjectMembershipDB]:
    """Get all members of a project with user details loaded."""
    stmt = (
        select(ProjectMembershipDB)
        .where(ProjectMembershipDB.project_id == project_id)
        .options(selectinload(ProjectMembershipDB.user))
    )
    result = await db.execute(stmt)
    return list(result.scalars().all())


def project_member_response_from_db(membership: ProjectMembershipDB) -> ProjectMemberResponse:
    """Helper function to convert ProjectMembershipDB to ProjectMemberResponse."""
    return ProjectMemberResponse(
        user_id=membership.user_id,
        user_name=membership.user.name,
        user_email=membership.user.email,
        user_is_active=membership.user.is_active,
        role=membership.role,
    )


async def get_project_membership(db: AsyncSession, project_id: int, user_id: int) -> ProjectMembershipDB:
    """Get a specific project membership by project ID and user ID."""
    stmt = select(ProjectMembershipDB).where(
        ProjectMembershipDB.project_id == project_id, ProjectMembershipDB.user_id == user_id
    )
    result = await db.execute(stmt)
    membership = result.scalar_one_or_none()
    if not membership:
        raise ProjectMemberNotFoundError()

    return membership


async def update_project_member_role(
    db: AsyncSession, project_id: int, user_id: int, new_role: ProjectRoles
) -> ProjectMembershipDB:
    """Update a project member's role."""
    membership = await get_project_membership(db, project_id, user_id)
    membership.role = new_role

    await db.commit()
    await db.refresh(membership)
    return membership


async def remove_project_member(db: AsyncSession, project_id: int, user_id: int) -> None:
    """Remove a user from a project."""
    membership = await get_project_membership(db, project_id, user_id)
    await db.delete(membership)
    await db.commit()


def has_required_role(user_role: ProjectRoles, required_role: ProjectRoles) -> bool:
    """
    Check if a user's role meets or exceeds the required role.
    (If you EDIT access, you also have READ access, etc.)
    """
    role_hierarchy = {
        ProjectRoles.READ: 1,
        ProjectRoles.EDIT: 2,
        ProjectRoles.MANAGE: 3,
    }
    return role_hierarchy[user_role] >= role_hierarchy[required_role]


async def _validate_project_create(db: AsyncSession, name: str, bucket_name: str) -> None:
    """
    Validate that both the project and bucket name are unique.
    """
    name_exists_stmt = select(exists().where(ProjectDB.name == name))
    name_exists_result = await db.execute(name_exists_stmt)
    if name_exists_result.scalar():
        raise ProjectCreationError(f"Project name '{name}' is already in use")

    bucket_exists_stmt = select(exists().where(ProjectDB.bucket_name == bucket_name))
    bucket_exists_result = await db.execute(bucket_exists_stmt)
    if bucket_exists_result.scalar():
        raise ProjectCreationError(f"Bucket name '{bucket_name}' is already in use")
