"""
Schemas for personal access token related endpoints.
"""

from pydantic import BaseModel, Field, model_validator
from typing_extensions import Self

# This needs to stay in sync with ProjectRoles in divbase-api.
# Redefined here to avoid a circular dependency on divbase-api.
ALLOWED_PROJECT_ROLES = ["read", "query", "edit", "manage"]


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

    @model_validator(mode="after")
    def no_project_scope_if_all_projects_true(self) -> Self:
        if self.all_projects and self.projects:
            raise ValueError("cannot specify project scope if 'all_projects' is True")
        return self

    @model_validator(mode="after")
    def validate_project_roles(self) -> Self:
        for project_id, role_str in self.projects.items():
            if role_str not in ALLOWED_PROJECT_ROLES:
                raise ValueError(f"invalid role '{role_str}' for project {project_id}")
        return self
