"""
CRUD operations for VCF dimensions for the Celery workers.
.
There are separate VCF dimensions CRUD functions for used with API endpoints in
packages/divbase-api/src/divbase_api/crud/vcf_dimensions.py
"""

import logging

from sqlalchemy import delete, select
from sqlalchemy.dialects.postgresql import insert
from sqlalchemy.orm.session import Session

from divbase_api.models.vcf_dimensions import SkippedVCFDB, VCFMetadataDB

logger = logging.getLogger(__name__)


def get_vcf_metadata_by_project(db: Session, project_id: int) -> dict:
    """
    FOR CELERY WORKERS, not for user interactions with API.

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
    FOR CELERY WORKERS, not for user interactions with API.

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
    FOR CELERY WORKERS, not for user interactions with API.

    Get VCF metadata by S3 key AND project ID (unique constraint).
    """
    stmt = select(VCFMetadataDB).where(
        VCFMetadataDB.vcf_file_s3_key == vcf_file_s3_key, VCFMetadataDB.project_id == project_id
    )
    result = db.execute(stmt)
    return result.scalar_one_or_none()


def create_or_update_vcf_metadata(db: Session, vcf_metadata_data: dict) -> None:
    """
    FOR CELERY WORKERS, not for user interactions with API.

    Upsert (insert or update) a VCF metadata entry in the database.

    This function uses PostgreSQL's ON CONFLICT DO UPDATE to ensure that if a VCF metadata entry
    with the same (vcf_file_s3_key, project_id) already exists, it will be updated with the new data.
    Otherwise, INSERT a new entry.
    """
    stmt = insert(VCFMetadataDB).values(**vcf_metadata_data)

    update_dict = {}
    for entry_index, entry_value in vcf_metadata_data.items():
        if entry_index not in ["vcf_file_s3_key", "project_id"]:
            update_dict[entry_index] = entry_value

    stmt = stmt.on_conflict_do_update(
        index_elements=["vcf_file_s3_key", "project_id"],
        set_=update_dict,
    )

    db.execute(stmt)
    db.commit()
    logger.info(
        f"VCF metadata created/updated for {vcf_metadata_data['vcf_file_s3_key']} "
        f"in project {vcf_metadata_data['project_id']}"
    )


def delete_vcf_metadata(db: Session, vcf_file_s3_key: str, project_id: int) -> None:
    """
    FOR CELERY WORKERS, not for user interactions with API.

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
    """
    FOR CELERY WORKERS, not for user interactions with API.

    Get skipped VCF entry by S3 key and project ID.
    """
    stmt = select(SkippedVCFDB).where(
        SkippedVCFDB.vcf_file_s3_key == vcf_file_s3_key, SkippedVCFDB.project_id == project_id
    )
    result = db.execute(stmt)
    return result.scalar_one_or_none()


def create_or_update_skipped_vcf(db: Session, skipped_vcf_data: dict) -> SkippedVCFDB:
    """
    FOR CELERY WORKERS, not for user interactions with API.

    Upsert (update or insert) skipped VCF entry. Similar to create_or_update_vcf_metadata but for tracking the skipped VCF files (=old divbase results VCF files).
    """
    stmt = insert(SkippedVCFDB).values(**skipped_vcf_data)

    update_dict = {}
    for entry_index, entry_value in skipped_vcf_data.items():
        if entry_index not in ["vcf_file_s3_key", "project_id"]:
            update_dict[entry_index] = entry_value

    stmt = stmt.on_conflict_do_update(
        index_elements=["vcf_file_s3_key", "project_id"],
        set_=update_dict,
    )

    db.execute(stmt)
    db.commit()

    logger.info(
        f"Skipped VCF entry created/updated for {skipped_vcf_data['vcf_file_s3_key']} in project {skipped_vcf_data['project_id']}"
    )


def delete_skipped_vcf(db: Session, vcf_file_s3_key: str, project_id: int) -> None:
    """
    FOR CELERY WORKERS, not for user interactions with API.

    Delete a skipped VCF entry by S3 key and project ID.
    """
    # TODO - add test for deleting a non existant entry, does it raise an error?
    stmt = delete(SkippedVCFDB).where(
        SkippedVCFDB.vcf_file_s3_key == vcf_file_s3_key, SkippedVCFDB.project_id == project_id
    )
    result = db.execute(stmt)
    db.commit()

    logger.info(
        f"Deleted skipped VCF entry for {vcf_file_s3_key} in project {project_id}. Rows affected: {result.rowcount}"
    )
