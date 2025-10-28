"""
API routes for VCF dimensions (technical metadata) operations.
"""

import logging

from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.vcf_dimensions import (
    create_or_update_vcf_metadata,
    get_vcf_metadata_by_keys,
    get_vcf_metadata_by_project,
)
from divbase_api.db import get_db
from divbase_api.deps import get_current_user
from divbase_api.models.users import UserDB

logger = logging.getLogger(__name__)

vcf_dimensions_router = APIRouter()


@vcf_dimensions_router.get("/projects/{project_id}")
async def list_vcf_metadata_for_project(
    project_id: int,
    db: AsyncSession = Depends(get_db),
    current_user: UserDB = Depends(get_current_user),
):
    """
    Get all VCF metadata entries for a project.
    Returns technical metadata (dimensions) for all VCF files in the project.
    """
    entries = await get_vcf_metadata_by_project(db, project_id)

    if not entries:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"VCF metadata not found for project {project_id}",
        )
    # TODO: is the detail giving away too much information?

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
                "indexed_at": entry.indexed_at.isoformat() if entry.indexed_at else None,
                "created_at": entry.created_at.isoformat() if entry.created_at else None,
                "updated_at": entry.updated_at.isoformat() if entry.updated_at else None,
            }
            for entry in entries
        ],
    }


@vcf_dimensions_router.get("/projects/{project_id}/files/{vcf_file_s3_key:path}")
async def get_vcf_metadata_for_specific_file_in_project(
    project_id: int,
    vcf_file_s3_key: str,
    db: AsyncSession = Depends(get_db),
    current_user: UserDB = Depends(get_current_user),
):
    """
    Get VCF metadata for a specific file in a project.
    Returns technical metadata (dimensions) for the specified VCF file.
    """
    entry = await get_vcf_metadata_by_keys(db, vcf_file_s3_key, project_id)

    if not entry:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"VCF metadata not found for file '{vcf_file_s3_key}' in project {project_id}",
        )
    # TODO: is the detail giving away too much information?

    return {
        "vcf_file_s3_key": entry.vcf_file_s3_key,
        "project_id": entry.project_id,
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


@vcf_dimensions_router.post("/", status_code=status.HTTP_201_CREATED)
async def create_or_update_vcf_metadata_entry(
    vcf_metadata_data: dict,
    db: AsyncSession = Depends(get_db),
    current_user: UserDB = Depends(get_current_user),
):
    """
    Create or update VCF metadata entry.

    This endpoint is typically called by worker tasks after indexing VCF files.
    Requires authentication.

    Request body should contain:
    - vcf_file_s3_key: str (required)
    - project_id: int (required)
    - s3_version_id: str
    - samples: list[str]
    - scaffolds: list[str]
    - variant_count: int
    - sample_count: int
    - file_size_bytes: int
    """

    if "vcf_file_s3_key" not in vcf_metadata_data or "project_id" not in vcf_metadata_data:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="vcf_file_s3_key and project_id are required",
        )

    try:
        entry = await create_or_update_vcf_metadata(db, vcf_metadata_data)

        return {
            "message": "VCF metadata created/updated successfully",
            "vcf_file_s3_key": entry.vcf_file_s3_key,
            "project_id": entry.project_id,
            "s3_version_id": entry.s3_version_id,
            "indexed_at": entry.indexed_at.isoformat() if entry.indexed_at else None,
            "created_at": entry.created_at.isoformat() if entry.created_at else None,
            "updated_at": entry.updated_at.isoformat() if entry.updated_at else None,
        }
    except Exception as e:
        logger.error(f"Error creating or inserting VCF metadata: {str(e)}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to create/update VCF metadata: {str(e)}",
        ) from e


# TODO add Delete entry route
