"""
Schemas for task history routes.
"""

from typing import Any, List, Optional

from pydantic import BaseModel, Field

# Request Models


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
    kwargs: Optional[str]
    eta: Optional[float]
    expires: Optional[float]
    retries: Optional[int]
    result: Optional[str]  # use raw string for now, but should convert to dict eventaully
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
