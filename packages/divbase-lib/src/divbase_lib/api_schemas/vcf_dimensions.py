"""
Schemas for VCF dimensions routes.
"""

from typing import Optional

from pydantic import BaseModel, Field


class DimensionUpdateKwargs(BaseModel):
    """Keyword arguments for dimension update task. Used to pass info to Celery task, and also for recording task history."""

    bucket_name: str
    project_id: int
    project_name: str
    user_id: int


class DimensionUpdateTaskResult(BaseModel):
    """Dimension update task result details. Based on the return of tasks.update_dimensions_index."""

    status: Optional[str] = None
    VCF_files_added: Optional[list[str]] = Field(
        None, description="VCF files that were added to dimensions index by this job"
    )
    VCF_files_skipped: Optional[list[str]] = Field(
        None, description="VCF files skipped by this job (previous DivBase-generated result VCFs)"
    )
    VCF_files_deleted: Optional[list[str]] = Field(
        None, description="VCF files that have been deleted from the project and thus have been dropped from the index"
    )


class DimensionsShowResult(BaseModel):
    """Result model for showing VCF dimensions for a project."""

    project_id: int
    project_name: str
    vcf_file_count: int
    vcf_files: list[dict]
    skipped_file_count: int
    skipped_files: list[dict]
