"""
Schemas for task history routes.
"""

from datetime import datetime
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
    created_at: datetime  # NOTE: This comes from TaskHistoryDB, so cannot be None unlike the other datetimes here.
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
    date_done: Optional[datetime] = None
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

    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    runtime: Optional[float] = None
