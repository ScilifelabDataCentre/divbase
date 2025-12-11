"""
Schemas for query routes.
"""

from typing import Any, Optional

from pydantic import BaseModel


# Request models
class SampleMetadataQueryRequest(BaseModel):
    """Request model for sample metadata query route."""

    tsv_filter: str
    metadata_tsv_name: str


class BcftoolsQueryRequest(BaseModel):
    """Request model for sample metadata query route."""

    tsv_filter: str
    metadata_tsv_name: str
    command: str  # TODO add field to decribe that this is bcftools commands


# Models for task kwargs and task results. Reused in task history schemas too, hence pydantic models and not just dataclasses.
class SampleMetadataQueryKwargs(BaseModel):
    """Keyword arguments for sample metadata query task. Used to pass info to Celery task, and also for recording task history."""

    tsv_filter: str
    metadata_tsv_name: str
    bucket_name: str
    project_id: int
    project_name: str
    user_id: int


class BcftoolsQueryKwargs(BaseModel):
    """Keyword arguments for BCFtools query task. Used to pass info to Celery task, and also for recording task history."""

    tsv_filter: str
    command: str
    metadata_tsv_name: str
    bucket_name: str
    project_id: int
    project_name: str
    user_id: int


class SampleMetadataQueryTaskResult(BaseModel):
    """Metadata query task result details. Based on the return of tasks.sample_metadata_query."""

    sample_and_filename_subset: list[dict[str, Any]]
    unique_sample_ids: list[str]
    unique_filenames: list[str]
    query_message: str
    status: Optional[str] = None


class BcftoolsQueryTaskResult(BaseModel):
    """BCFtools query task result details. Based on the return of tasks.bcftools_query."""

    output_file: str
    status: Optional[str] = None
