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
            "description": (
                "At most one of tsv_filter or samples may be provided. "
                "If neither is provided, all project samples/files are used."
            )
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

    @model_validator(mode="after")
    def validate_sample_selection_mode(self) -> Self:
        tsv_filter = self.tsv_filter
        samples = self.samples

        if tsv_filter is not None and samples is not None:
            raise ValueError("Only one of tsv_filter or samples may be provided.")

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

    @model_validator(mode="after")
    def validate_sample_selection_mode(self) -> Self:
        tsv_filter = self.tsv_filter
        samples = self.samples

        if tsv_filter is not None and samples is not None:
            raise ValueError("Only one of tsv_filter or samples may be provided.")

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
