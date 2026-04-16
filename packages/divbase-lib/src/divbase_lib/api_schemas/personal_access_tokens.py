"""
Schemas for personal access token related endpoints.
"""

from pydantic import BaseModel, Field


class PATPermissions(BaseModel):
    """
    Personal Access Token permissions model.
    This is stored in the db as JSONB in the PersonalAccessTokenDB.permissions.
    When null, the PAT has same access as underlying user.

    Parsed from the JSONB stored in the PAT DB column via model_validate.
    If the DB column is NULL, the PAT has full, unrestricted access (equivalent to a JWT token).

    - all_projects: True = access all the user's projects with same role as user
                  False = access only the projects listed in `projects`.
    - projects: {str(project_id): role_string}; only enforced when all_projects=False.
    - task_history: endpoint level scope;. For the 2 task history routes not tied to a specific project.
        (these are task_history for a specific task id and for the user)
    """

    all_projects: bool = False
    projects: dict[str, str] = Field(default_factory=dict)
    task_history: bool = False
