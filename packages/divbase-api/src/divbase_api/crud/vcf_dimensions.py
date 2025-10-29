"""
CRUD operations for VCF dimensions.
"""

import logging

from sqlalchemy import delete, select
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.models.vcf_dimensions import SkippedVCFDB, VCFMetadataDB

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


async def delete_vcf_metadata(db: AsyncSession, vcf_file_s3_key: str, project_id: int) -> None:
    """
    Delete a VCF metadata entry by S3 key and project ID.

    Called when a VCF file is removed from the project bucket.
    """
    stmt = delete(VCFMetadataDB).where(
        VCFMetadataDB.vcf_file_s3_key == vcf_file_s3_key, VCFMetadataDB.project_id == project_id
    )
    await db.execute(stmt)
    await db.commit()

    logger.info(f"Deleted VCF metadata for {vcf_file_s3_key} in project {project_id}.")


async def get_skipped_vcf_by_keys(db: AsyncSession, vcf_file_s3_key: str, project_id: int) -> SkippedVCFDB | None:
    """Get skipped VCF entry by S3 key and project ID."""
    stmt = select(SkippedVCFDB).where(
        SkippedVCFDB.vcf_file_s3_key == vcf_file_s3_key, SkippedVCFDB.project_id == project_id
    )
    result = await db.execute(stmt)
    return result.scalar_one_or_none()


async def create_or_update_skipped_vcf(db: AsyncSession, skipped_vcf_data: dict) -> SkippedVCFDB:
    """
    Upsert (update or insert) skipped VCF entry.

    Called by workers when a VCF file is identified as a DivBase result file based on its header.

    Args:
        db: Database session
        skipped_vcf_data: Dictionary with fields:
            - vcf_file_s3_key: str (required)
            - project_id: int (required)
            - s3_version_id: str
            - skip_reason: str (e.g., "divbase_generated")
    """
    existing_entry = await get_skipped_vcf_by_keys(
        db, skipped_vcf_data["vcf_file_s3_key"], skipped_vcf_data["project_id"]
    )

    if existing_entry:
        for key, value in skipped_vcf_data.items():
            if key not in ["vcf_file_s3_key", "project_id"]:
                setattr(existing_entry, key, value)
        skipped_entry = existing_entry
        logger.info(f"Updated skipped VCF entry for {skipped_vcf_data['vcf_file_s3_key']}")
    else:
        skipped_entry = SkippedVCFDB(**skipped_vcf_data)
        db.add(skipped_entry)
        logger.info(f"Created skipped VCF entry for {skipped_vcf_data['vcf_file_s3_key']}")

    await db.commit()
    await db.refresh(skipped_entry)
    return skipped_entry


async def get_skipped_vcfs_by_project(db: AsyncSession, project_id: int) -> list[SkippedVCFDB]:
    """Get all skipped VCF entries for a given project."""
    stmt = select(SkippedVCFDB).where(SkippedVCFDB.project_id == project_id)
    result = await db.execute(stmt)
    return list(result.scalars().all())


async def delete_skipped_vcf(db: AsyncSession, vcf_file_s3_key: str, project_id: int) -> None:
    """Delete a skipped VCF entry by S3 key and project ID."""
    stmt = delete(SkippedVCFDB).where(
        SkippedVCFDB.vcf_file_s3_key == vcf_file_s3_key, SkippedVCFDB.project_id == project_id
    )
    result = await db.execute(stmt)
    await db.commit()

    logger.info(
        f"Deleted skipped VCF entry for {vcf_file_s3_key} in project {project_id}. Rows affected: {result.rowcount}"
    )
