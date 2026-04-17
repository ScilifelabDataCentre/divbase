"""
Schemas for personal access token related endpoints.
"""

from pydantic import BaseModel, Field


class PATPermissions(BaseModel):
    """
    Personal Access Token permissions model.
    This is stored in the db as JSONB in the PersonalAccessTokenDB.permissions.

    - all_projects: True = access all the user's projects with same role as user
                    False = access only the projects listed in `projects`.
    - projects: {str(project_id): role_string}; only enforced when all_projects=False.
        Where role_string must be a valid ProjectRoles enum value.
    - task_history: endpoint level scope.
        For the 2 task history routes not tied to a specific project.
        (these are 'task-history id' and for 'task-history user')
    """

    all_projects: bool = False
    projects: dict[str, str] = Field(default_factory=dict)
    task_history: bool = False

    # TODO - model validate, can't have projects and all_projects True at the same time,
    # and role strings must be valid ProjectRoles enum values.
