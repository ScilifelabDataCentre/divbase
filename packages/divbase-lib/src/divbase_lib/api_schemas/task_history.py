"""
Schemas for task history routes.
"""

from typing import Any, Optional, Union

from pydantic import BaseModel, Field

from divbase_lib.api_schemas.queries import (
    BcftoolsQueryKwargs,
    BcftoolsQueryTaskResult,
    SampleMetadataQueryKwargs,
    SampleMetadataQueryTaskResult,
)
from divbase_lib.api_schemas.vcf_dimensions import DimensionUpdateKwargs, DimensionUpdateTaskResult

# Response Models


class FlowerTaskResult(BaseModel):
    """Task details as returned by the Flower API"""

    uuid: str
    name: Optional[str] = None
    state: Optional[str] = None
    received: Optional[float] = None
    sent: Optional[float] = None
    started: Optional[float] = None
    rejected: Optional[float] = None
    succeeded: Optional[float] = None
    failed: Optional[float] = None
    retried: Optional[float] = None
    revoked: Optional[float] = None
    args: Optional[str] = None
    kwargs: Optional[Union[SampleMetadataQueryKwargs, BcftoolsQueryKwargs, DimensionUpdateKwargs, dict[str, Any]]] = (
        None
    )
    eta: Optional[float] = None
    expires: Optional[float] = None
    retries: Optional[int] = None
    result: Optional[Union[BcftoolsQueryTaskResult, SampleMetadataQueryTaskResult, DimensionUpdateTaskResult, str]] = (
        None
    )
    exception: Optional[str] = None
    timestamp: Optional[float] = None
    runtime: Optional[float] = None
    traceback: Optional[str] = None
    exchange: Optional[str] = None
    routing_key: Optional[str] = None
    clock: Optional[int] = None
    client: Optional[str] = None
    root: Optional[str] = None
    root_id: Optional[str] = None
    parent: Optional[str] = None
    parent_id: Optional[str] = None
    children: Optional[list[Any]] = None
    worker: Optional[str] = None


class TaskHistoryResults(BaseModel):
    """Results from all tasks fetched from a task history request."""

    tasks: dict[str, FlowerTaskResult] = Field(..., description="Mapping of task_id to Flower task details")
    user_email: Optional[str] = Field(None, description="Email of the user who requested the task history")
