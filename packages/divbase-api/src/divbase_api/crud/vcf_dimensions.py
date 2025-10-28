"""
CRUD operations for VCF dimensions.
"""

import logging

from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.models.vcf_dimensions import VCFMetadataDB

logger = logging.getLogger(__name__)


async def get_vcf_metadata_by_keys(db: AsyncSession, vcf_file_s3_key: str, project_id: int) -> VCFMetadataDB | None:
    """
    Get VCF metadata by S3 key AND project ID (unique constraint).
    """
    stmt = select(VCFMetadataDB).where(
        VCFMetadataDB.vcf_file_s3_key == vcf_file_s3_key, VCFMetadataDB.project_id == project_id
    )
    result = await db.execute(stmt)
    return result.scalar_one_or_none()


async def create_or_update_vcf_metadata(db: AsyncSession, vcf_metadata_data: dict) -> VCFMetadataDB:
    """
    Upsert (update or insert) VCF metadata entry.
    Called by workers after processing VCF file dimensions.

    Args:
        db: Database session
        vcf_metadata_data: Dictionary with VCF metadata fields:
            - vcf_file_s3_key: str (required)
            - project_id: int (required)
            - s3_version_id: str
            - samples: list[str]
            - scaffolds: list[str]
            - variant_count: int
            - sample_count: int
            - file_size_bytes: int
    """
    existing_entry = await get_vcf_metadata_by_keys(
        db, vcf_metadata_data["vcf_file_s3_key"], vcf_metadata_data["project_id"]
    )

    if existing_entry:
        for key, value in vcf_metadata_data.items():
            if key not in ["vcf_file_s3_key", "project_id"]:
                setattr(existing_entry, key, value)
        dimensions_entry = existing_entry
        logger.info(f"Updated VCF metadata for {vcf_metadata_data['vcf_file_s3_key']}")
    else:
        dimensions_entry = VCFMetadataDB(**vcf_metadata_data)
        db.add(dimensions_entry)
        logger.info(f"Created VCF metadata for {vcf_metadata_data['vcf_file_s3_key']}")

    await db.commit()
    await db.refresh(dimensions_entry)
    return dimensions_entry


async def get_vcf_metadata_by_project(db: AsyncSession, project_id: int) -> list[VCFMetadataDB]:
    """
    Get all VCF metadata entries for a given project ID.
    """
    stmt = select(VCFMetadataDB).where(VCFMetadataDB.project_id == project_id)
    result = await db.execute(stmt)
    return list(result.scalars().all())


# TODO: Delete entry when file is removed

# TODO: Delete all entries for a project

# TODO: get vcf metadata entries by project ID (not vcf_file_s3_key)

# TODO: get skipped VCF results file by keys - should that perhaps be its own table?
