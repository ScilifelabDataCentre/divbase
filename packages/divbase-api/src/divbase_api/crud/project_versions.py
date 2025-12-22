"""
CRUD operations for project versions.

Project versions correspond to the overall state of all files in a project's bucket
at a given timestamp. Each version entry records the mapping of file names to their unique
version IDs (hashes) at that point in time.

Version entries are created and managed via the API.
"""

import logging
from datetime import datetime, timezone

from fastapi.concurrency import run_in_threadpool
from sqlalchemy import select
from sqlalchemy.exc import IntegrityError
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.exceptions import (
    ProjectVersionAlreadyExistsError,
    ProjectVersionCreationError,
    ProjectVersionNotFoundError,
)
from divbase_api.models.project_versions import ProjectVersionDB
from divbase_api.models.projects import ProjectDB
from divbase_api.services.s3_client import S3FileManager
from divbase_lib.api_schemas.project_versions import (
    AddVersionResponse,
    DeleteVersionResponse,
    ProjectVersionDetailResponse,
    ProjectVersionInfo,
)

logger = logging.getLogger(__name__)


async def add_project_version(
    db: AsyncSession,
    name: str,
    description: str | None,
    project: ProjectDB,
    user_id: int,
    s3_file_manager: S3FileManager,
) -> AddVersionResponse:
    """Add a new project version entry into the database."""
    files = await run_in_threadpool(s3_file_manager.latest_version_of_all_files, bucket_name=project.bucket_name)

    if not files:
        raise ProjectVersionCreationError(
            message="Cannot create version entry for an empty project. No files uploaded to the project yet.",
        )

    new_version = ProjectVersionDB(
        name=name,
        description=description,
        files=files,
        project_id=project.id,
        user_id=user_id,
    )

    db.add(new_version)
    try:
        await db.commit()
    except IntegrityError as e:
        await db.rollback()
        error_details = str(e.orig).lower()
        if "unique_name_project" in error_details and "unique constraint" in error_details:
            raise ProjectVersionAlreadyExistsError(
                message=f"A project version with the name '{name}' already exists for this project. Please choose a different name.",
            ) from None
        else:
            logger.error(f"Unexpected db integrity error when attempting to add a new project version: {e}")
            raise

    await db.refresh(new_version)
    return AddVersionResponse(
        name=new_version.name, description=new_version.description, created_at=new_version.created_at.isoformat()
    )


async def list_project_versions(
    db: AsyncSession,
    project_id: int,
    include_deleted: bool = False,
) -> list[ProjectVersionInfo]:
    """List all version entries for a given project."""
    stmt = select(
        ProjectVersionDB.name,
        ProjectVersionDB.description,
        ProjectVersionDB.created_at,
        ProjectVersionDB.is_deleted,
    )
    stmt = stmt.where(ProjectVersionDB.project_id == project_id)
    if not include_deleted:
        stmt = stmt.where(ProjectVersionDB.is_deleted == False)  # noqa: E712
    stmt = stmt.order_by(ProjectVersionDB.created_at.desc())

    result = await db.execute(stmt)
    versions = [
        ProjectVersionInfo(
            name=row.name,
            description=row.description,
            created_at=row.created_at.isoformat(),
            is_deleted=row.is_deleted,
        )
        for row in result
    ]
    return versions


async def get_project_version_details(
    db: AsyncSession,
    project_id: int,
    version_name: str,
) -> ProjectVersionDetailResponse:
    """Get the mapping of files to their version IDs at a specific project version."""
    stmt = select(
        ProjectVersionDB.name,
        ProjectVersionDB.description,
        ProjectVersionDB.created_at,
        ProjectVersionDB.is_deleted,
        ProjectVersionDB.files,
    )
    stmt = stmt.where(
        ProjectVersionDB.project_id == project_id,
        ProjectVersionDB.name == version_name,
    )

    result = await db.execute(stmt)
    version_entry = result.one_or_none()
    if version_entry is None:
        raise ProjectVersionNotFoundError(message=f"Version '{version_name}' not found for the project.")

    name, description, created_at, is_deleted, files = version_entry
    return ProjectVersionDetailResponse(
        name=name,
        description=description,
        created_at=created_at.isoformat(),
        is_deleted=is_deleted,
        files=files,
    )


async def soft_delete_version(
    db: AsyncSession,
    project_id: int,
    version_name: str,
) -> DeleteVersionResponse:
    """
    Soft delete a project version entry.
    A cron job is responsible for permanently deleting soft deleted versions after a grace period.
    """
    stmt = select(ProjectVersionDB).where(
        ProjectVersionDB.project_id == project_id,
        ProjectVersionDB.name == version_name,
    )
    result = await db.execute(stmt)
    version_entry = result.scalar_one_or_none()
    if version_entry is None:
        raise ProjectVersionNotFoundError(
            message=f"Version '{version_name}' not found for the project. Perhaps it was already deleted or you mistyped the version name?"
        )

    if version_entry.is_deleted and version_entry.date_deleted is not None:
        return DeleteVersionResponse(
            name=version_name,
            already_deleted=True,
            date_deleted=version_entry.date_deleted.isoformat(),
        )

    version_entry.is_deleted = True
    version_entry.date_deleted = datetime.now(tz=timezone.utc)
    await db.commit()
    await db.refresh(version_entry)
    return DeleteVersionResponse(
        name=version_name,
        already_deleted=False,
        date_deleted=version_entry.date_deleted.isoformat(),
    )
