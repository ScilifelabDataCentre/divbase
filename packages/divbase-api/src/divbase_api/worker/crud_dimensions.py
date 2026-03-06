"""
CRUD operations for VCF dimensions for the Celery workers
There are separate VCF dimensions CRUD functions for used with API endpoints in
packages/divbase-api/src/divbase_api/crud/vcf_dimensions.py
"""

import dataclasses
import logging
from dataclasses import dataclass
from typing import List

from sqlalchemy import delete, select
from sqlalchemy.dialects.postgresql import insert
from sqlalchemy.orm import selectinload
from sqlalchemy.orm.session import Session

from divbase_api.models.vcf_dimensions import SkippedVCFDB, VCFMetadataDB, VCFMetadataSamplesDB, VCFMetadataScaffoldsDB

logger = logging.getLogger(__name__)


@dataclass
class VCFMetadataData:
    vcf_file_s3_key: str
    project_id: int
    s3_version_id: str | None
    samples: List[str]
    scaffolds: List[str]
    variant_count: int
    sample_count: int
    file_size_bytes: int


@dataclass
class SkippedVCFData:
    vcf_file_s3_key: str
    project_id: int
    s3_version_id: str | None
    skip_reason: str


def get_vcf_metadata_by_project(db: Session, project_id: int) -> dict:
    """
    FOR CELERY WORKERS, not for user interactions with API.

    Get all VCF metadata entries for a given project ID.
    """
    stmt = (
        select(VCFMetadataDB)
        .where(VCFMetadataDB.project_id == project_id)
        .options(selectinload(VCFMetadataDB.samples), selectinload(VCFMetadataDB.scaffolds))
    )
    db_result = db.execute(stmt)

    entries = list(db_result.scalars().all())

    result = {
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


def create_or_update_vcf_metadata(db: Session, vcf_metadata_data: VCFMetadataData) -> None:
    """
    FOR CELERY WORKERS, not for user interactions with API.

    Upsert (insert or update) a VCF metadata entry in the database.

    This function uses PostgreSQL's ON CONFLICT DO UPDATE to ensure that if a VCF metadata entry
    with the same (vcf_file_s3_key, project_id) already exists, it will be updated with the new data.
    Otherwise, INSERT a new entry.

    Samples and scaffolds are stored in separate FK tables (VCFMetadataSamplesDB, VCFMetadataScaffoldsDB)
    and are popped (cut out) from the dict before the upsert to VCFMetadataDB and then added to the FK tables accordingly.
    """
    data_dict = dataclasses.asdict(vcf_metadata_data)
    samples: list[str] = data_dict.pop("samples", [])
    scaffolds: list[str] = data_dict.pop("scaffolds", [])

    stmt = insert(VCFMetadataDB).values(**data_dict)

    update_dict = {k: v for k, v in data_dict.items() if k not in ["vcf_file_s3_key", "project_id"]}

    stmt = stmt.on_conflict_do_update(
        index_elements=["vcf_file_s3_key", "project_id"],
        set_=update_dict,
    )

    db.execute(stmt)

    # Flush instead of commit here so that the INSERT/UPDATE is visible to get_vcf_metadata_by_keys below.
    # The flush will be visible to other operations within the same db session, but will not yet be committed to the database (invisivle to other db sessions).
    # This avoids concurrency issues for the child tables if two identical dimensions update jobs are run concurrently in the job system.
    db.flush()

    # Get the upserted/updated entry from the main model (=parent object) and use it to add the samples and scaffolds to the respective FK tables (=child objects)
    vcf_metadata = get_vcf_metadata_by_keys(db, vcf_metadata_data.vcf_file_s3_key, vcf_metadata_data.project_id)
    # Important! This is a full replacement based on the VCF files in the bucket, not an append. If a samples in a vcf file is added, changed, or removed, the relationship will delete the existing samples entries for that VCF with cascade="all, delete-orphan" and insert the new ones.
    vcf_metadata.samples = [VCFMetadataSamplesDB(sample_name=name) for name in samples]
    vcf_metadata.scaffolds = [VCFMetadataScaffoldsDB(scaffold_name=name) for name in scaffolds]

    # Single commit for parent upsert and child table updates.
    db.commit()

    logger.info(
        f"VCF Dimensions created/updated for {vcf_metadata_data.vcf_file_s3_key} "
        f"in project {vcf_metadata_data.project_id}"
    )


def delete_vcf_metadata(db: Session, vcf_file_s3_key: str, project_id: int) -> None:
    """
    FOR CELERY WORKERS, not for user interactions with API.

    Delete a VCF metadata entry by S3 key and project ID. The batch version in delete_vcf_metadata_batch
    is more efficient and is preferred over this. This function is kept in the codebase for legacy and edge case reason.

    Called when a VCF file is removed from the project bucket.
    """
    stmt = delete(VCFMetadataDB).where(
        VCFMetadataDB.vcf_file_s3_key == vcf_file_s3_key, VCFMetadataDB.project_id == project_id
    )
    result = db.execute(stmt)
    db.commit()

    logger.info(f"Deleted VCF metadata for {vcf_file_s3_key} in project {project_id}. Rows affected: {result.rowcount}")


def delete_vcf_metadata_batch(db: Session, vcf_file_s3_key_batch: list[str], project_id: int) -> None:
    """
    FOR CELERY WORKERS, not for user interactions with API.

    Batch delete VCF metadata entries by a list of S3 keys and project ID.
    This is a version of delete_vcf_metadata that delete multiple entries in a single transaction.

    Called when multiple VCF files are removed from the project bucket, e.g. when an entire project is deleted.
    """
    stmt = delete(VCFMetadataDB).where(
        VCFMetadataDB.vcf_file_s3_key.in_(vcf_file_s3_key_batch), VCFMetadataDB.project_id == project_id
    )
    result = db.execute(stmt)
    db.commit()

    logger.info(
        f"Batch deleted VCF metadata for {len(vcf_file_s3_key_batch)} files in project {project_id}. "
        f"Rows affected: {result.rowcount}"
    )


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


def create_or_update_skipped_vcf(db: Session, skipped_vcf_data: SkippedVCFData) -> None:
    """
    FOR CELERY WORKERS, not for user interactions with API.

    Upsert (update or insert) skipped VCF entry. Similar to create_or_update_vcf_metadata but for tracking the skipped VCF files (=old divbase results VCF files).
    """
    data_dict = dataclasses.asdict(skipped_vcf_data)

    stmt = insert(SkippedVCFDB).values(**data_dict)

    update_dict = {}
    for entry_index, entry_value in data_dict.items():
        if entry_index not in ["vcf_file_s3_key", "project_id"]:
            update_dict[entry_index] = entry_value

    stmt = stmt.on_conflict_do_update(
        index_elements=["vcf_file_s3_key", "project_id"],
        set_=update_dict,
    )

    db.execute(stmt)
    db.commit()

    logger.info(
        f"Skipped VCF entry created/updated for {skipped_vcf_data.vcf_file_s3_key} in project {skipped_vcf_data.project_id}"
    )


def delete_skipped_vcf(db: Session, vcf_file_s3_key: str, project_id: int) -> None:
    """
    FOR CELERY WORKERS, not for user interactions with API.

    Delete a skipped VCF entry by S3 key and project ID. The batch version in delete_skipped_vcf_batch
    is more efficient and is preferred over this. This function is kept in the codebase for legacy and edge case reasons.
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


def delete_skipped_vcf_batch(db: Session, vcf_file_s3_key_batch: list[str], project_id: int) -> None:
    """
    FOR CELERY WORKERS, not for user interactions with API.

    Batch delete skipped VCF entries by a list of S3 keys and project ID.
    This is a version of delete_skipped_vcf that deletes multiple entries in a single transaction.

    Called when multiple skipped VCF files are removed from the project bucket.
    """
    stmt = delete(SkippedVCFDB).where(
        SkippedVCFDB.vcf_file_s3_key.in_(vcf_file_s3_key_batch), SkippedVCFDB.project_id == project_id
    )
    result = db.execute(stmt)
    db.commit()

    logger.info(
        f"Batch deleted skipped VCF entries for {len(vcf_file_s3_key_batch)} files in project {project_id}. "
        f"Rows affected: {result.rowcount}"
    )
