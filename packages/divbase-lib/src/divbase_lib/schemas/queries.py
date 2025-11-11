"""
Schemas for query routes.
"""

from typing import Any, Dict, List

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
    user_name: str


class BcftoolsQueryKwargs(BaseModel):
    """Keyword arguments for BCFtools query task. Used to pass info to Celery task, and also for recording task history."""

    tsv_filter: str
    command: str
    metadata_tsv_name: str
    bucket_name: str
    project_id: int
    user_name: str


class TaskMetadataQueryResult(BaseModel):
    """Metadata query task result details. Based on the return of tasks.sample_metadata_query."""

    sample_and_filename_subset: List[Dict[str, Any]]
    unique_sample_ids: List[str]
    unique_filenames: List[str]
    query_message: str
    status: str


# TODO these task results models could be defined elsewhere and used here and by the routes for queries, dimensions, and task_hist
class TaskBcftoolsQueryResult(BaseModel):
    """BCFtools query task result details. Based on the return of tasks.bcftools_query."""

    status: str
    output_file: str
    submitter: str
