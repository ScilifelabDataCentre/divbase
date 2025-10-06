"""
S3 files routes. Provides a mix of direct responses and pre-signed URLs depending on the operation.
"""

import logging

from fastapi import APIRouter, Query, status

from divbase_api.schemas.files import DownloadObjectResponse, DownloadObjectsRequest
from divbase_api.services.pre_signed_urls import S3PreSignedService

logger = logging.getLogger(__name__)

files_router = APIRouter()


@files_router.get("/download", status_code=status.HTTP_200_OK, response_model=list[DownloadObjectResponse])
async def download_files(
    bucket_name: str = "local-project-1",
    object_names: list[str] = Query(..., description="List of object names to download"),
):
    """Download a file from S3."""
    s3_signer_service = S3PreSignedService()

    response = []
    for obj_name in object_names:
        url = s3_signer_service.create_presigned_url_for_download(
            bucket_name=bucket_name, object_name=obj_name, version_id=None
        )
        response.append(
            DownloadObjectResponse(
                object_name=obj_name,
                pre_signed_url=url,
                version_id=None,
            )
        )
    return response


# Post request instead of GET as get doesn't support/encourage body content.
@files_router.post("/download_at_version", status_code=status.HTTP_200_OK, response_model=list[DownloadObjectResponse])
async def download_files_at_version(
    objects_list: DownloadObjectsRequest,
    bucket_name: str = "local-project-1",
):
    """Download a files at specific versions from S3."""
    s3_signer_service = S3PreSignedService()

    response = []
    for obj in objects_list.objects:
        url = s3_signer_service.create_presigned_url_for_download(
            bucket_name=bucket_name, object_name=obj.object_name, version_id=obj.version_id
        )
        response.append(
            DownloadObjectResponse(
                object_name=obj.object_name,
                pre_signed_url=url,
                version_id=obj.version_id,
            )
        )

    return response
