import logging

from sqlalchemy import delete, select
from sqlalchemy.orm.session import Session

from divbase_api.models.vcf_dimensions import SkippedVCFDB, VCFMetadataDB

logger = logging.getLogger(__name__)


def get_vcf_metadata_by_project(db: Session, project_id: int) -> dict:
    """
    Get all VCF metadata entries for a given project ID.
    """
    stmt = select(VCFMetadataDB).where(VCFMetadataDB.project_id == project_id)
    db_result = db.execute(stmt)

    entries = list(db_result.scalars().all())

    result = {
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

    return result


def get_skipped_vcfs_by_project_worker(db: Session, project_id: int) -> dict[str, str]:
    """
    Get all skipped VCF entries for a given project.
    """
    stmt = select(SkippedVCFDB).where(SkippedVCFDB.project_id == project_id)
    result = db.execute(stmt)
    entries = list(result.scalars().all())

    already_skipped_vcfs = {}
    for entry in entries:
        already_skipped_vcfs[entry.vcf_file_s3_key] = entry.s3_version_id

    return already_skipped_vcfs


def get_vcf_metadata_by_keys(db: Session, vcf_file_s3_key: str, project_id: int) -> VCFMetadataDB | None:
    """
    Get VCF metadata by S3 key AND project ID (unique constraint).
    """
    stmt = select(VCFMetadataDB).where(
        VCFMetadataDB.vcf_file_s3_key == vcf_file_s3_key, VCFMetadataDB.project_id == project_id
    )
    result = db.execute(stmt)
    return result.scalar_one_or_none()


def create_or_update_vcf_metadata(db: Session, vcf_metadata_data: dict) -> None:
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
    existing_entry = get_vcf_metadata_by_keys(db, vcf_metadata_data["vcf_file_s3_key"], vcf_metadata_data["project_id"])

    if existing_entry:
        for key, value in vcf_metadata_data.items():
            if key not in ["vcf_file_s3_key", "project_id"]:
                setattr(existing_entry, key, value)
        dimensions_entry = existing_entry
    else:
        dimensions_entry = VCFMetadataDB(**vcf_metadata_data)
        db.add(dimensions_entry)

    db.commit()
    db.refresh(dimensions_entry)
    logger.info(
        f"VCF metadata created/updated for {dimensions_entry.vcf_file_s3_key} in project {dimensions_entry.project_id}"
    )


def delete_vcf_metadata(db: Session, vcf_file_s3_key: str, project_id: int) -> None:
    """
    Delete a VCF metadata entry by S3 key and project ID.

    Called when a VCF file is removed from the project bucket.
    """
    stmt = delete(VCFMetadataDB).where(
        VCFMetadataDB.vcf_file_s3_key == vcf_file_s3_key, VCFMetadataDB.project_id == project_id
    )
    result = db.execute(stmt)
    db.commit()

    logger.info(f"Deleted VCF metadata for {vcf_file_s3_key} in project {project_id}. Rows affected: {result.rowcount}")


def get_skipped_vcf_by_keys(db: Session, vcf_file_s3_key: str, project_id: int) -> SkippedVCFDB | None:
    """Get skipped VCF entry by S3 key and project ID."""
    stmt = select(SkippedVCFDB).where(
        SkippedVCFDB.vcf_file_s3_key == vcf_file_s3_key, SkippedVCFDB.project_id == project_id
    )
    result = db.execute(stmt)
    return result.scalar_one_or_none()


def create_or_update_skipped_vcf(db: Session, skipped_vcf_data: dict) -> SkippedVCFDB:
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
    existing_entry = get_skipped_vcf_by_keys(db, skipped_vcf_data["vcf_file_s3_key"], skipped_vcf_data["project_id"])

    if existing_entry:
        for key, value in skipped_vcf_data.items():
            if key not in ["vcf_file_s3_key", "project_id"]:
                setattr(existing_entry, key, value)
        skipped_entry = existing_entry
    else:
        skipped_entry = SkippedVCFDB(**skipped_vcf_data)
        db.add(skipped_entry)

    db.commit()
    db.refresh(skipped_entry)
    logger.info(
        f"Skipped VCF entry {'created' if not existing_entry else 'updated'} for {skipped_entry.vcf_file_s3_key} in project {skipped_entry.project_id}"
    )
    return skipped_entry


def delete_skipped_vcf(db: Session, vcf_file_s3_key: str, project_id: int) -> None:
    """Delete a skipped VCF entry by S3 key and project ID."""
    # TODO - add test for deleting a non existant entry, does it raise an error?
    stmt = delete(SkippedVCFDB).where(
        SkippedVCFDB.vcf_file_s3_key == vcf_file_s3_key, SkippedVCFDB.project_id == project_id
    )
    result = db.execute(stmt)
    db.commit()

    logger.info(
        f"Deleted skipped VCF entry for {vcf_file_s3_key} in project {project_id}. Rows affected: {result.rowcount}"
    )
