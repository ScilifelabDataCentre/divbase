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
from divbase_lib.api_schemas.vcf_dimensions import DimensionUpdateTaskResult


class TaskHistoryResult(BaseModel):
    """
    Task details as returned by queries to the SQAlchemy+pg results backend.
    """

    uuid: str
    submitter_email: Optional[str] = None
    status: Optional[str] = None
    result: Optional[
        Union[
            dict[
                str, Any
            ],  # Note! This dict must come first here so that error results are preserved and not incorrectly inserted into the result models
            SampleMetadataQueryTaskResult,
            BcftoolsQueryTaskResult,
            DimensionUpdateTaskResult,
        ]
    ] = None
    date_done: Optional[str] = None
    name: Optional[str] = None
    args: Optional[str] = None
    kwargs: Optional[
        Union[
            SampleMetadataQueryKwargs,
            BcftoolsQueryKwargs,
        ]
    ] = None
    worker: Optional[str] = None
    created_at: Optional[float] = None
    started_at: Optional[float] = None
    completed_at: Optional[float] = None
    runtime: Optional[float] = None


# TODO consider if traceback: Optional[str] = None;     retries: Optional[int] = None;     queue: Optional[str] = None would be relevant?


class TaskHistoryResults(BaseModel):
    """Results from all tasks fetched from a task history request."""

    tasks: dict[str, TaskHistoryResult] = Field(..., description="Mapping of task_id to results backend task details")
    user_email: Optional[str] = Field(None, description="Email of the user who requested the task history")
