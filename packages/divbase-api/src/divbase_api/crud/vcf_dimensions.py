"""
CRUD operations for VCF dimensions.

The functions in this file are intended to be used with API endpoints.
There are separate VCF dimensions CRUD functions for the Celery workers in
packages/divbase-api/src/divbase_api/worker/crud_dimensions.py
"""

import structlog
from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.orm import selectinload

from divbase_api.exceptions import DimensionsUpdateAlreadyInProcessError
from divbase_api.models.task_history import CeleryTaskMeta, TaskHistoryDB
from divbase_api.models.vcf_dimensions import SkippedVCFDB, VCFMetadataDB, VCFMetadataSamplesDB, VCFMetadataScaffoldsDB

logger = structlog.get_logger(__name__)


async def get_vcf_metadata_by_project_async(db: AsyncSession, project_id: int) -> dict:
    """
    FOR USER INTERACTIONS WITH API ONLY
    Get all VCF metadata entries for a given project ID, including the related sample and scaffold names.

    VCFMetadataDB has a one-to-many relationship relationship with VCFMetadataSamplesDB and with VCFMetadataScaffoldsDB.

    Eager load with selectinload is used to minimize the number of db queries (to 3 queries). It is used instead of table joins since these child tables can have a large amount of entries
    per parent VCF file, which would result in a large amount of duplicated data being loaded if table joins were used and be inefficient.
    In this eager load, 1 query is used to load all VCFMetadataDB entries (=VCF files) for the project, 1 query is used to load all related VCFMetadataSamplesDB entries (sample names)
    for those VCFMetadataDB entries, and 1 query is used to load all related VCFMetadataScaffoldsDB entries (scaffold names) for those VCFMetadataDB entries.

    Lazy load, which is not used here, would instead result in 1 + 2N queries for N VCFMetadataDB entries (=VCF files) for the project, since each VCFMetadataDB entry would be queried
    separately for its related sample and scaffold names.
    """

    stmt = (
        select(VCFMetadataDB)
        .where(VCFMetadataDB.project_id == project_id)
        .options(selectinload(VCFMetadataDB.samples), selectinload(VCFMetadataDB.scaffolds))
    )
    result = await db.execute(stmt)
    entries = list(result.scalars().all())

    return {
        "project_id": project_id,
        "vcf_file_count": len(entries),
        "vcf_files": [
            {
                "vcf_file_s3_key": entry.vcf_file_s3_key,
                "s3_version_id": entry.s3_version_id,
                "samples": [s.sample_name for s in entry.samples],
                "scaffolds": [s.scaffold_name for s in entry.scaffolds],
                "variant_count": entry.variant_count,
                "sample_count": entry.sample_count,
                "file_size_bytes": entry.file_size_bytes,
                "created_at": entry.created_at.isoformat(),
                "updated_at": entry.updated_at.isoformat(),
            }
            for entry in entries
        ],
    }


async def get_skipped_vcfs_by_project_async(db: AsyncSession, project_id: int) -> list[SkippedVCFDB]:
    """
    FOR USER INTERACTIONS WITH API ONLY, Celery workers have their own dimensions CRUD functions.

    Get all skipped VCF entries for a given project.
    """
    stmt = select(SkippedVCFDB).where(SkippedVCFDB.project_id == project_id)
    result = await db.execute(stmt)
    return list(result.scalars().all())


async def get_unique_samples_by_project_async(db: AsyncSession, project_id: int) -> list[str]:
    """
    Get unique sample names across all VCF files from a project's dimensions entries.
    """

    stmt = (
        select(VCFMetadataSamplesDB.sample_name)
        .join(VCFMetadataDB, VCFMetadataSamplesDB.vcf_metadata_id == VCFMetadataDB.id)
        .where(VCFMetadataDB.project_id == project_id)
        .distinct()
        .order_by(VCFMetadataSamplesDB.sample_name)
    )
    result = await db.execute(stmt)
    return list(result.scalars().all())


async def get_unique_scaffolds_by_project_async(db: AsyncSession, project_id: int) -> list[str]:
    """
    Get unique scaffold names across all VCF files for a project.
    """

    stmt = (
        select(VCFMetadataScaffoldsDB.scaffold_name)
        .join(VCFMetadataDB, VCFMetadataScaffoldsDB.vcf_metadata_id == VCFMetadataDB.id)
        .where(VCFMetadataDB.project_id == project_id)
        .distinct()
        .order_by(VCFMetadataScaffoldsDB.scaffold_name)
    )
    result = await db.execute(stmt)
    scaffolds = list(result.scalars().all())

    # Sort scaffold names in the same way as the dimensions show CLI does when returning all dimensions data: numeric first, then alphabetic
    # Numeric sorting of name strings results means that 10 comes after 2 for scaffolds that have numeric names.
    numeric = sorted([int(s) for s in scaffolds if s.isdigit()])
    non_numeric = sorted([s for s in scaffolds if not s.isdigit()])
    return [str(n) for n in numeric] + non_numeric


async def get_unique_vcf_files_by_project_async(db: AsyncSession, project_id: int) -> list[dict[str, str]]:
    """
    Get unique VCF file/version pairs across all VCF files for a project.
    """

    stmt = (
        select(VCFMetadataDB.vcf_file_s3_key, VCFMetadataDB.s3_version_id)
        .where(VCFMetadataDB.project_id == project_id)
        .distinct()
        .order_by(VCFMetadataDB.vcf_file_s3_key)
    )
    result = await db.execute(stmt)
    rows = result.all()
    unique_vcf_files = []
    for file_name, version_id in rows:
        unique_vcf_files.append({"vcf_file_s3_key": file_name, "s3_version_id": version_id})
    return unique_vcf_files


async def check_no_dimensions_update_task_already_in_progress(
    db: AsyncSession, project_id: int, project_name: str
) -> None:
    """
    Check that there is not already a dimensions update job in process for the given project and raise if so.
    This check includes jobs that are both queued and running.

    There only needs to be one dimensions update job run at a time per project, otherwise it is just waited compute + bandwidth.
    """
    ongoing_dimensions_tasks_subq = (
        select(CeleryTaskMeta.task_id)
        .where(CeleryTaskMeta.status.in_(["PENDING", "STARTED"]))
        .where(CeleryTaskMeta.name == "tasks.update_vcf_dimensions_task")
        .subquery()
    )

    # the .id is the user facing id (aka rolling int) and not the celery internal uuid which is .task_id
    stmt = (
        select(TaskHistoryDB.id)
        .join(ongoing_dimensions_tasks_subq, TaskHistoryDB.task_id == ongoing_dimensions_tasks_subq.c.task_id)
        .where(TaskHistoryDB.project_id == project_id)
    )
    result = await db.execute(stmt)
    # NOTE: don't use one_or_none() here since we can't guarantee there is only one ongoing task.
    ongoing_task_id = result.first()

    if ongoing_task_id is not None:
        raise DimensionsUpdateAlreadyInProcessError(
            project_name=project_name,
            ongoing_task_id=ongoing_task_id[0],
        )
