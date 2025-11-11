"""
Schemas for task history routes.
"""

from typing import Any, Dict, List, Optional, Union

from pydantic import BaseModel, Field

# Response Models


class TaskMetadataQueryResult(BaseModel):
    """Metadata query task result details. Based on the return of tasks.sample_metadata_query."""

    sample_and_filename_subset: List[Dict[str, Any]]
    unique_sample_ids: List[str]
    unique_filenames: List[str]
    query_message: str
    status: str


class TaskBcftoolsQueryResult(BaseModel):
    """BCFtools query task result details. Based on the return of tasks.bcftools_query."""

    status: str
    output_file: str
    submitter: str


class TaskDimensionUpdateResult(BaseModel):
    """Dimension update task result details. Based on the return of tasks.update_dimensions_index."""

    status: str
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


class FlowerTaskResult(BaseModel):
    """Task details as returned by the Flower API"""

    uuid: str
    name: Optional[str]
    state: Optional[str]
    received: Optional[float]
    sent: Optional[float]
    started: Optional[float]
    rejected: Optional[float]
    succeeded: Optional[float]
    failed: Optional[float]
    retried: Optional[float]
    revoked: Optional[float]
    args: Optional[str]
    kwargs: Optional[str]
    eta: Optional[float]
    expires: Optional[float]
    retries: Optional[int]
    result: Optional[Union[TaskBcftoolsQueryResult, TaskMetadataQueryResult, TaskDimensionUpdateResult, str]]
    exception: Optional[str]
    timestamp: Optional[float]
    runtime: Optional[float]
    traceback: Optional[str]
    exchange: Optional[str]
    routing_key: Optional[str]
    clock: Optional[int]
    client: Optional[str]
    root: Optional[str]
    root_id: Optional[str]
    parent: Optional[str]
    parent_id: Optional[str]
    children: Optional[List[Any]]
    worker: Optional[str]


class TaskHistoryResults(BaseModel):
    """Results from all tasks fetched from a task history request."""

    tasks: dict[str, FlowerTaskResult] = Field(..., description="Mapping of task_id to Flower task details")
    user_email: Optional[str] = Field(None, description="Email of the user who requested the task history")
