"""
Schemas for task history routes.
"""

from typing import Any, Dict, List, Optional, Union

from pydantic import BaseModel, Field

from divbase_lib.schemas.queries import (
    BcftoolsQueryKwargs,
    BcftoolsQueryTaskResult,
    SampleMetadataQueryKwargs,
    SampleMetadataQueryTaskResult,
)
from divbase_lib.schemas.vcf_dimensions import DimensionUpdateKwargs, DimensionUpdateTaskResult

# Response Models


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
    kwargs: Optional[Union[SampleMetadataQueryKwargs, BcftoolsQueryKwargs, DimensionUpdateKwargs, Dict[str, Any]]]
    eta: Optional[float]
    expires: Optional[float]
    retries: Optional[int]
    result: Optional[Union[BcftoolsQueryTaskResult, SampleMetadataQueryTaskResult, DimensionUpdateTaskResult, str]]
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
