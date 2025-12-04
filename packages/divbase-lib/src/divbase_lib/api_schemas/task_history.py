"""
Schemas for task history routes.
"""

from datetime import datetime
from typing import Any, Optional, Union

from pydantic import BaseModel, Field

from divbase_lib.api_schemas.queries import (
    BcftoolsQueryKwargs,
    BcftoolsQueryTaskResult,
    SampleMetadataQueryKwargs,
    SampleMetadataQueryTaskResult,
)
from divbase_lib.api_schemas.vcf_dimensions import DimensionUpdateKwargs, DimensionUpdateTaskResult


class TaskHistoryResult(BaseModel):
    """
    Task details as returned by queries to the SQAlchemy+pg results backend.
    """

    uuid: str
    user_id: Optional[int] = None
    project_id: Optional[int] = None
    submitter_email: str
    status: Optional[str] = None
    created_at: datetime
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    runtime: Optional[float] = None
    result: Optional[
        Union[
            BcftoolsQueryTaskResult,
            SampleMetadataQueryTaskResult,
            DimensionUpdateTaskResult,
            dict[str, Any],  # For error results
            str,
        ]
    ] = None
    date_done: Optional[datetime] = None
    name: Optional[str] = None
    args: Optional[str] = None
    kwargs: Optional[Union[SampleMetadataQueryKwargs, BcftoolsQueryKwargs, DimensionUpdateKwargs, dict[str, Any]]] = (
        None
    )
    worker: Optional[str] = None


# TODO consider if traceback: Optional[str] = None;     retries: Optional[int] = None;     queue: Optional[str] = None would be relevant?


class TaskHistoryResults(BaseModel):
    """Results from all tasks fetched from a task history request."""

    tasks: dict[str, TaskHistoryResult] = Field(..., description="Mapping of task_id to results backend task details")
    user_email: Optional[str] = Field(None, description="Email of the user who requested the task history")
