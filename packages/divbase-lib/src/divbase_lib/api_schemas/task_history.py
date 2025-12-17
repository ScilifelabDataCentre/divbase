"""
Schemas for task history routes.
"""

from typing import Any, Optional, Union

from pydantic import BaseModel

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

    id: int
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
            DimensionUpdateKwargs,
        ]
    ] = None
    worker: Optional[str] = None
    created_at: Optional[str] = None
    started_at: Optional[str] = None
    completed_at: Optional[str] = None
    runtime: Optional[float] = None
