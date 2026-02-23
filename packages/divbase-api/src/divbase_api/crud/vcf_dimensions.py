"""
CRUD operations for VCF dimensions.
"""

import logging

from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.orm import selectinload

from divbase_api.models.vcf_dimensions import SkippedVCFDB, VCFMetadataDB, VCFMetadataSamplesDB, VCFMetadataScaffoldsDB

logger = logging.getLogger(__name__)


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
    FOR USER INTERACTIONS WITH API ONLY
    Get all skipped VCF entries for a given project.
    """
    stmt = select(SkippedVCFDB).where(SkippedVCFDB.project_id == project_id)
    result = await db.execute(stmt)
    return list(result.scalars().all())


async def get_unique_samples_by_project_async(db: AsyncSession, project_id: int) -> list[str]:
    """
    Get unique sample names across all VCF files from a project's dimensions entries.

    Samples are stored in as ARRAY(String) in the VCFMetadataDB model and need to be flattened before finding the unqiue values.
    To do all operations on the PostgreSQL side (to avoid having do it here in the fastAPI side), need to first use unnest() to flatten the arrays.Expand commentComment on lines R58 to R59Resolved
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

    Like samples, scaffolds are stored in as ARRAY(String) in the VCFMetadataDB model and need to be flattened with the unnest() PostgreSQL function.
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
