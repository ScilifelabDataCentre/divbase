"""
Schemas for task history routes.
"""

from typing import Any, Dict, List, Optional, Union

from pydantic import BaseModel, Field

# Response Models


class FlowerTaskMetadataQueryResult(BaseModel):
    sample_and_filename_subset: Optional[list[Dict[str, Any]]]
    unique_sample_ids: Optional[List[str]]
    unique_filenames: Optional[List[str]]
    query_message: Optional[str]
    status: Optional[str]
    task_id: Optional[str]


class FlowerTaskBcftoolsQueryResult(BaseModel):
    status: str
    output_file: Optional[str]
    submitter: Optional[str]


class FlowerTaskDimensionUpdateResult(BaseModel):
    status: str
    submitter: Optional[str]
    VCF_files_added: Optional[List[str]] = Field(
        None, alias="VCF files that were added to dimensions index by this job"
    )
    VCF_files_skipped: Optional[List[str]] = Field(
        None, alias="VCF files skipped by this job (previous DivBase-generated result VCFs)"
    )
    VCF_files_deleted: Optional[List[str]] = Field(
        None, alias="VCF files that have been deleted from the project and thus have been dropped from the index"
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
    result: Optional[
        Union[FlowerTaskBcftoolsQueryResult, FlowerTaskMetadataQueryResult, FlowerTaskDimensionUpdateResult, str]
    ]
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
