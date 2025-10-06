"""
S3 files routes. Provides a mix of direct responses and pre-signed URLs depending on the operation.
"""

import logging

from fastapi import APIRouter, Query, status

from divbase_api.config import settings
from divbase_api.schemas.files import DownloadObjectResponse, DownloadObjectsRequest
from divbase_api.services.pre_signed_urls import S3PreSignedService
from divbase_lib.s3_client import S3FileManager

logger = logging.getLogger(__name__)

files_router = APIRouter()


@files_router.get("/download", status_code=status.HTTP_200_OK, response_model=list[DownloadObjectResponse])
async def download_files(
    project_name: str = "local-project-1",
    object_names: list[str] = Query(..., description="List of object names to download"),
):
    """Download a file from S3."""
    bucket_name = project_name  # TODO, in future handle mapping

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
    project_name: str = "local-project-1",
):
    """Download a files at specific versions from S3."""
    bucket_name = project_name  # TODO, in future handle mapping
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


@files_router.get("/list", status_code=status.HTTP_200_OK)
async def list_files(
    project_name: str = "local-project-1",
):
    """List all files in the project's bucket"""
    bucket_name = project_name  # TODO, in future handle mapping

    s3_file_manager = S3FileManager(
        url=settings.s3.endpoint_url,
        access_key=settings.s3.access_key.get_secret_value(),
        secret_key=settings.s3.secret_key.get_secret_value(),
    )

    return s3_file_manager.list_files(bucket_name=bucket_name)
