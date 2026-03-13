"""
Single source-of-truth for Celery task names. Stored in separate file to avoid circular imports
"""

from enum import Enum


class TaskName(str, Enum):
    """
    Single source-of-truth for Celery task names. Allows these names name to be reused in other layer of the codebase.
    """

    SAMPLE_METADATA_QUERY = "tasks.sample_metadata_query"
    BCFTOOLS_QUERY = "tasks.bcftools_query"
    UPDATE_VCF_DIMENSIONS = "tasks.update_vcf_dimensions_task"
