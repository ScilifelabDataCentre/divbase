import logging

from sqlalchemy import select
from sqlalchemy.orm.session import Session

from divbase_api.models.vcf_dimensions import VCFMetadataDB

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
                "indexed_at": entry.indexed_at.isoformat() if entry.indexed_at else None,
                "created_at": entry.created_at.isoformat() if entry.created_at else None,
                "updated_at": entry.updated_at.isoformat() if entry.updated_at else None,
            }
            for entry in entries
        ],
    }

    if not result["vcf_files"]:
        raise FileNotFoundError("No VCF metadata found for project {project_id}")

    return result
