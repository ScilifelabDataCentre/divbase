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
    duplicate_of_job_id: Optional[int] = Field(
        None,
        description="The job ID of the task that this task is a duplicate of. Only set if the task is a duplicate.",
    )
    message: Optional[str] = Field(
        None, description="A message describing the outcome of the task. Only set if the task is a duplicate."
    )


class DimensionsUpdateSubmitResult(BaseModel):
    """Result model for submitting a dimensions update job."""

    job_id: int
    outcome: str = Field(
        ...,
        description="Whether a new job was enqueued or an existing active job was reused. Allowed values: 'new', 'existing'.",
    )


class DimensionsShowResult(BaseModel):
    """Result model for showing VCF dimensions for a project."""

    project_id: int
    project_name: str
    vcf_file_count: int
    vcf_files: list[dict]
    skipped_file_count: int
    skipped_files: list[dict]


class DimensionsSamplesResult(BaseModel):
    """Result model for showing unique samples across project VCFs."""

    unique_samples: list[str]  # Already sorted, by the CRUD function get_unique_samples_by_project_async()


class DimensionsScaffoldsResult(BaseModel):
    """Result model for showing unique scaffolds across project VCFs."""

    unique_scaffolds: list[str]  # Already sorted, by the CRUD function get_unique_scaffolds_by_project_async()
