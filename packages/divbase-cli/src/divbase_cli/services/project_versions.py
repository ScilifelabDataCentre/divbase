"""
Service layer for DivBase CLI project version operations.
"""

from divbase_cli.user_auth import make_authenticated_request
from divbase_lib.api_schemas.project_versions import (
    AddVersionRequest,
    AddVersionResponse,
    DeleteVersionRequest,
    DeleteVersionResponse,
    ProjectVersionDetailResponse,
    ProjectVersionInfo,
)


def add_version_command(project_name: str, divbase_base_url: str, name: str, description: str) -> AddVersionResponse:
    """Add a new version to the project versions table stored on the divbase server"""
    request_data = AddVersionRequest(name=name, description=description)

    response = make_authenticated_request(
        method="PATCH",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/project-versions/add?project_name={project_name}",
        json=request_data.model_dump(),
    )

    return AddVersionResponse(**response.json())


def list_versions_command(project_name: str, include_deleted: bool, divbase_base_url: str) -> list[ProjectVersionInfo]:
    """
    List all versions in the project versions table stored on the divbase server.
    Returns a dict of version names (keys) to details about the versions.
    """
    response = make_authenticated_request(
        method="GET",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/project-versions/list?project_name={project_name}&include_deleted={str(include_deleted).lower()}",
    )

    project_versions = []
    response_data = response.json()
    for version in response_data:
        project_versions.append(ProjectVersionInfo(**version))

    return project_versions


def get_version_details_command(
    project_name: str, divbase_base_url: str, version_name: str
) -> ProjectVersionDetailResponse:
    """Get details about a specific project version, including all files and their version IDs at that version."""
    response = make_authenticated_request(
        method="GET",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/project-versions/version_details?project_name={project_name}&version_name={version_name}",
    )

    return ProjectVersionDetailResponse(**response.json())


def delete_version_command(project_name: str, divbase_base_url: str, version_name: str) -> DeleteVersionResponse:
    """
    Delete a version from the project versions table stored on the divbase server.
    This marks the version as (soft) deleted server side,
    and it will eventually be permanently deleted (after some grace period).
    """
    request_data = DeleteVersionRequest(version_name=version_name)

    response = make_authenticated_request(
        method="DELETE",
        divbase_base_url=divbase_base_url,
        api_route=f"v1/project-versions/delete?project_name={project_name}",
        json=request_data.model_dump(),
    )

    return DeleteVersionResponse(**response.json())
