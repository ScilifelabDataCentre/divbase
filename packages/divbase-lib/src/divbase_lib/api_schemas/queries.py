"""
Schemas for query routes.
"""

from typing import Optional, Self

from pydantic import BaseModel, ConfigDict, field_validator, model_validator

from divbase_lib.utils import split_semicolon_bcftools_command_segments

# Shared helper functions and models


def validate_command_not_empty(command: str) -> str:
    """
    Validator function to ensure that the user-submitted bcftools command string is not empty or just whitespace.
    If there are multiple segments separated by semicolons, it also validates that none of the segments are empty or whitespace.

    Shared helper for command field validators in BcftoolsQueryRequest.command and BcftoolsQueryKwargs.command.
    """
    command_string = command.strip()
    if command_string == "":
        raise ValueError(
            "The --command option must be a non-empty bcftools view string. "
            'If you only want to subset based on samples, use --command "view -s" in combination with one of the sample-selection options: --tsv-filter, --samples, --samples-file.'
        )

    segments = split_semicolon_bcftools_command_segments(command_string)
    for position, segment in enumerate(segments, start=1):
        if segment.strip() == "":
            raise ValueError(
                f"The --command option has an empty pipeline segment at position {position}. "
                "Use semicolons only between complete bcftools view commands."
            )

    return command


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
    command: str  # bcftools command input with --command
    samples: list[str] | None = None
    all_samples: bool = False

    @field_validator("command")
    @classmethod
    def validate_command(_cls, value: str) -> str:
        return validate_command_not_empty(value)

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
    command: str  # bcftools command input with --command
    bucket_name: str
    project_id: int
    project_name: str
    user_id: int
    job_id: int
    samples: list[str] | None = None
    all_samples: bool = False

    @field_validator("command")
    @classmethod
    def validate_command(_cls, value: str) -> str:
        return validate_command_not_empty(value)

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
            # To ensure that task-history deserialization does not break for existing tasks.
            # Could be handled by a backfilling migration instead.
            self.all_samples = True
            return self

        if tsv_filter is not None and self.metadata_tsv_name is None:
            raise ValueError("metadata_tsv_name must be provided when tsv_filter is used.")

        if samples is not None and len(samples) == 0:
            raise ValueError("samples must contain at least one sample ID when provided.")

        return self


class SampleFileMappingResult(BaseModel):
    sample_id: str
    filename: str

    @model_validator(mode="before")
    @classmethod
    def map_legacy_task_history_keys(cls, data):
        """
        Support historical task-history payloads that stored uppercase keys.
        """
        if not isinstance(data, dict):
            return data

        normalized = dict(data)
        if "sample_id" not in normalized:
            if "Sample_ID" in normalized:
                normalized["sample_id"] = normalized["Sample_ID"]
            elif "#Sample_ID" in normalized:
                normalized["sample_id"] = normalized["#Sample_ID"]
        if "filename" not in normalized and "Filename" in normalized:
            normalized["filename"] = normalized["Filename"]

        return normalized


class SampleMetadataQueryTaskResult(BaseModel):
    """Metadata query task result details. Based on the return of tasks.sample_metadata_query."""

    sample_and_filename_subset: list[SampleFileMappingResult]
    unique_sample_ids: list[str]
    unique_filenames: list[str]
    query_message: str
    warnings: list[str] = []
    status: str | None = None


class BcftoolsQueryTaskResult(BaseModel):
    """BCFtools query task result details. Based on the return of tasks.bcftools_query."""

    output_file: str
    status: Optional[str] = None
