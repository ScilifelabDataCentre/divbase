"""
Schemas for VCF dimensions routes.
"""

from typing import List, Optional

from pydantic import BaseModel, Field


class DimensionUpdateKwargs(BaseModel):
    """Keyword arguments for dimension update task. Used to pass info to Celery task, and also for recording task history."""

    bucket_name: str
    project_id: int
    user_name: str


class DimensionUpdateTaskResult(BaseModel):
    """Dimension update task result details. Based on the return of tasks.update_dimensions_index."""

    status: Optional[str] = None
    submitter: str
    VCF_files_added: Optional[List[str]] = Field(
        None, description="VCF files that were added to dimensions index by this job"
    )
    VCF_files_skipped: Optional[List[str]] = Field(
        None, description="VCF files skipped by this job (previous DivBase-generated result VCFs)"
    )
    VCF_files_deleted: Optional[List[str]] = Field(
        None, description="VCF files that have been deleted from the project and thus have been dropped from the index"
    )


class DimensionsShowResult(BaseModel):
    """Result model for showing VCF dimensions for a project."""

    project_id: int
    project_name: str
    vcf_file_count: int
    vcf_files: List[dict]
    skipped_file_count: int
    skipped_files: List[dict]
