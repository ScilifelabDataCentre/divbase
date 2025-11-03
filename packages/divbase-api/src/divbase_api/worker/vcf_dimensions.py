import logging

from fastapi import APIRouter
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.vcf_dimensions import (
    get_vcf_metadata_by_project,
)

logger = logging.getLogger(__name__)

vcf_dimensions_router = APIRouter()


async def list_vcf_metadata_for_project(
    project_id: int,
    db: AsyncSession,
):
    """
    Get all VCF metadata entries for a project.

    Requires worker service account authentication.
    """
    result = await get_vcf_metadata_by_project(db, project_id)

    if not result["vcf_files"]:
        raise FileNotFoundError("No VCF metadata found for project {project_id}")

    return result
