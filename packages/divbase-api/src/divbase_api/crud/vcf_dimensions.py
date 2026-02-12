"""
CRUD operations for VCF dimensions.
"""

import logging

from sqlalchemy import func, select
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.models.vcf_dimensions import SkippedVCFDB, VCFMetadataDB

logger = logging.getLogger(__name__)


async def get_vcf_metadata_by_project_async(db: AsyncSession, project_id: int) -> dict:
    """
    FOR USER INTERACTIONS WITH API ONLY
    Get all VCF metadata entries for a given project ID.
    """
    stmt = select(VCFMetadataDB).where(VCFMetadataDB.project_id == project_id)
    result = await db.execute(stmt)
    entries = list(result.scalars().all())

    return {
        "project_id": project_id,
        "vcf_file_count": len(entries),
        "vcf_files": [
            {
                "vcf_file_s3_key": entry.vcf_file_s3_key,
                "s3_version_id": entry.s3_version_id,
                "samples": entry.samples,
                "scaffolds": entry.scaffolds,
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
    To do all operations on the PostgreSQL side (to avoid having do it here in the fastAPI side), need to first use unnest() to flatten the arrays.
    """

    stmt = select(func.unnest(VCFMetadataDB.samples)).where(VCFMetadataDB.project_id == project_id).distinct()
    result = await db.execute(stmt)
    return sorted([row[0] for row in result])
