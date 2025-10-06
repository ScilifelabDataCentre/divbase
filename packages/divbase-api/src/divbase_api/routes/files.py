"""
S3 files routes. Provides a mix of direct responses and pre-signed URLs depending on the operation.
"""

import logging

from fastapi import APIRouter, HTTPException, Query, status

from divbase_api.services.pre_signed_urls import S3PreSignedService

logger = logging.getLogger(__name__)

files_router = APIRouter()


@files_router.get("/download", status_code=status.HTTP_200_OK)
async def download_files(
    bucket_name: str = "local-project-1",
    object_names: list[str] = Query(..., description="List of object names to download"),
):
    """Download a file from S3."""
    s3_signer_service = S3PreSignedService()
    dload_urls = s3_signer_service.make_download_urls(bucket_name=bucket_name, object_names=object_names)

    if not any(dload_urls.values()):
        raise HTTPException(
            status_code=404, detail="No file(s) found, check the bucket and/or object names are correct."
        )
    missing_files = [name for name, url in dload_urls.items() if not url]
    return {"pre_signed_urls": dload_urls, "missing_files": missing_files if missing_files else None}
