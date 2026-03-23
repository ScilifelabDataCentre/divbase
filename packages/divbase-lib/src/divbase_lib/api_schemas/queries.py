"""
Schemas for query routes.
"""

from typing import Any, Optional, Self

from pydantic import BaseModel, ConfigDict, model_validator


class SharedBaseModel(BaseModel):
    """Shared pydantic BaseModel for VCF query response and kwarg schemas."""

    model_config = ConfigDict(
        validate_assignment=True,
        json_schema_extra={
            "description": ("Exactly one sample-selection mode must be provided: tsv_filter, samples, or all_samples.")
        },
    )


# Request models
class SampleMetadataQueryRequest(BaseModel):
    """Request model for sample metadata query route."""

    tsv_filter: str
    metadata_tsv_name: str


class BcftoolsQueryRequest(SharedBaseModel):
    """
    Request model for sample metadata query route.
    See ConfigDict in SharedBaseModel for validation rules around tsv_filter and samples fields.
    """

    tsv_filter: str | None = None  # Used for metadata mode only
    metadata_tsv_name: str | None = None  # Used for metadata mode only
    command: str  # TODO add field to describe that this is bcftools commands
    samples: list[str] | None = None
    all_samples: bool = False

    @model_validator(mode="after")
    def validate_sample_selection_mode(self) -> Self:
        tsv_filter = self.tsv_filter
        samples = self.samples
        all_samples = self.all_samples

        selection_count = 0
        if tsv_filter is not None:
            selection_count += 1
        if samples is not None:
            selection_count += 1
        if all_samples:
            selection_count += 1

        if selection_count > 1:
            raise ValueError("Only one of tsv_filter, samples, or all_samples may be provided.")
        if selection_count == 0:
            raise ValueError("One sample-selection mode must be provided (tsv_filter, samples, or all_samples).")

        if tsv_filter is not None and self.metadata_tsv_name is None:
            raise ValueError("metadata_tsv_name must be provided when tsv_filter is used.")

        if samples is not None and len(samples) == 0:
            raise ValueError("samples must contain at least one sample ID when provided.")

        return self


# Models for task kwargs and task results. Reused in task history schemas too, hence pydantic models and not just dataclasses.
class SampleMetadataQueryKwargs(BaseModel):
    """Keyword arguments for sample metadata query task. Used to pass info to Celery task, and also for recording task history."""

    tsv_filter: str
    metadata_tsv_name: str
    bucket_name: str
    project_id: int
    project_name: str
    user_id: int


class BcftoolsQueryKwargs(SharedBaseModel):
    """
    Keyword arguments for BCFtools query task. Used to pass info to Celery task, and also for recording task history.
    See ConfigDict in SharedBaseModel for validation rules around tsv_filter and samples fields.
    """

    tsv_filter: str | None = None  # Used for metadata mode only
    metadata_tsv_name: str | None = None  # Used for metadata mode only
    command: str
    bucket_name: str
    project_id: int
    project_name: str
    user_id: int
    job_id: int
    samples: list[str] | None = None
    all_samples: bool = False

    @model_validator(mode="after")
    def validate_sample_selection_mode(self) -> Self:
        tsv_filter = self.tsv_filter
        samples = self.samples
        all_samples = self.all_samples

        selection_count = 0
        if tsv_filter is not None:
            selection_count += 1
        if samples is not None:
            selection_count += 1
        if all_samples:
            selection_count += 1

        if selection_count > 1:
            raise ValueError("Only one of tsv_filter, samples, or all_samples may be provided.")
        if selection_count == 0:
            # Backward compatibility for historical task kwargs created before all_samples option was implemented.
            # To ensure that task-history deserizliation does not break for existing tasks.
            # Could be handled by a backfilling migration instead.
            self.all_samples = True
            return self

        if tsv_filter is not None and self.metadata_tsv_name is None:
            raise ValueError("metadata_tsv_name must be provided when tsv_filter is used.")

        if samples is not None and len(samples) == 0:
            raise ValueError("samples must contain at least one sample ID when provided.")

        return self


class SampleMetadataQueryTaskResult(BaseModel):
    """Metadata query task result details. Based on the return of tasks.sample_metadata_query."""

    sample_and_filename_subset: list[dict[str, Any]]
    unique_sample_ids: list[str]
    unique_filenames: list[str]
    query_message: str
    warnings: list[str] = []
    status: str | None = None


class BcftoolsQueryTaskResult(BaseModel):
    """BCFtools query task result details. Based on the return of tasks.bcftools_query."""

    output_file: str
    status: Optional[str] = None
