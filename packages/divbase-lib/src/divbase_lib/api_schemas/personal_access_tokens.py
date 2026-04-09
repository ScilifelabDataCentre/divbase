"""
Schemas for personal access token related endpoints.
"""

from datetime import datetime

from pydantic import AwareDatetime, BaseModel, ConfigDict, Field


class PATPermissions(BaseModel):
    """
    Personal Access Token permissions model.
    This is stored in the db as JSONB in the PersonalAccessTokenDB.permissions.
    Field is nullable — when null, the PAT has same access as underlying user.

    Parsed from the JSONB stored in the PAT DB column via model_validate.
    If the DB column is NULL, the PAT has full, unrestricted access (equivalent to a JWT token).

    - all_projects: True = access all the user's projects with same role as user
                  False = access only the projects listed in `projects`.
    - projects: {str(project_id): role_string}; only enforced when all_projects=False.
    - task_history / whoami: endpoint/route level scopes as not tied to a specific project.
    """

    all_projects: bool = True
    projects: dict[str, str] = Field(default_factory=dict)
    task_history: bool = False
    whoami: bool = False


class PATCreateRequest(BaseModel):
    name: str = Field(..., min_length=1, max_length=100, description="Name of the PAT")
    description: str | None = Field(None, max_length=500, description="Optional description of the PAT.")
    permissions: PATPermissions | None = Field(
        None,
        description="Optional scoped permissions of the PAT, If None the PAT has the same access as the user.",
    )
    expires_at: AwareDatetime | None = Field(None, description="Optional expiration date for the PAT")


class PATCreateResponse(BaseModel):
    id: int = Field(..., description="ID of the PAT.")
    name: str = Field(..., description="Name of the PAT.")
    description: str | None = Field(None, description="Description of the PAT.")
    permissions: PATPermissions | None = Field(None, description="Scoped permissions for this PAT.")
    expires_at: AwareDatetime | None = Field(None, description="Expiration date for the PAT.")
    token: str = Field(..., description="The PAT in plaintext. Only shown at creation time. Store it securely.")

    model_config = ConfigDict(from_attributes=True)


class PATItem(BaseModel):
    """A single PAT item in a list of PATs returned to a user."""

    id: int
    name: str = Field(..., description="Name of the PAT, for user reference.")
    description: str | None = Field(None, description="Optional description of the PAT.")
    permissions: PATPermissions | None = Field(None, description="Scoped permissions for this PAT.")
    expires_at: AwareDatetime | None = Field(None, description="Expiration date for the PAT.")
    last_used_at: datetime | None = Field(None, description="Last used date for the PAT.")
    created_at: datetime = Field(..., description="Creation date for the PAT.")

    model_config = ConfigDict(from_attributes=True)
