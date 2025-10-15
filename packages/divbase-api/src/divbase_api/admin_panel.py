"""
Manage the starlette-admin interface for the DivBase API.

The views created for each model rely on overriding some of the default behavior provided by starlette-admin.
These overrides are on methods inside BaseModelView (parent of ModelView, which each custom class is inheriting from).

"""

from typing import Any

from fastapi import FastAPI
from pydantic import SecretStr
from sqlalchemy.ext.asyncio import AsyncEngine
from starlette.requests import Request
from starlette_admin import BooleanField, EmailField, IntegerField, StringField, TextAreaField
from starlette_admin.contrib.sqla import Admin, ModelView

from divbase_api.models.projects import ProjectDB, ProjectMembershipDB
from divbase_api.models.users import UserDB
from divbase_api.security import get_password_hash

PAGINATION_DEFAULTS = [5, 10, 25, -1]  # (for number of items per page toggle)


class UserView(ModelView):
    """
    Custom admin panel View for the UserDB model.

    Some Notes:
    - To handle password hashing, we override the create method to hash the password before calling the original create method.
    - We disable deletion of users, you can mark a user as soft deleted though.
    - Adding/editing project memberships is disabled in this view, they should be handled in the ProjectMembership view.
    - Changing password for existing user is not supported: We will support a custom password reset flow later.
    """

    page_size_options = PAGINATION_DEFAULTS
    fields = [
        "id",
        StringField("name", required=True),
        EmailField("email", required=True),
        StringField("password", required=True),
        StringField(
            "hashed_password",
            required=False,
            disabled=True,
            help_text="Hashed password is auto created by the system from password, cannot be edited.",
        ),
        "is_admin",
        "is_active",
        "is_deleted",
        "project_memberships",
    ]

    exclude_fields_from_list = ["hashed_password", "password"]
    exclude_fields_from_create = [
        "project_memberships",
        "is_deleted",
        "is_active",
    ]  # hashed_password wont be defined yet but needs to be included.
    exclude_fields_from_edit = ["project_memberships", "password", "hashed_password"]
    exclude_fields_from_detail = ["hashed_password", "password"]

    def can_delete(self, request: Request) -> bool:
        return False

    async def create(self, request: Request, data: dict) -> Any:
        """
        We override the default create method so we can take a password from the frontend form
        and create a hashed password which will be needed to add to the db.
        """
        password = data["password"]
        hashed_password = get_password_hash(SecretStr(password))
        data["hashed_password"] = hashed_password

        return await super().create(request, data)


class ProjectView(ModelView):
    """
    Custom admin panel View for the ProjectDB model.

    Project memberships are managed in the ProjectMembership view.
    """

    page_size_options = PAGINATION_DEFAULTS
    fields = [
        "id",
        StringField("name", required=True, label="Project Name"),
        TextAreaField("description", required=False, label="Description"),
        StringField("bucket_name", required=True, label="Bucket Name"),
        IntegerField(
            "storage_quota_bytes",
            required=True,
            label="Storage Quota (Bytes)",
            help_text="Maximum storage allowed for this project in bytes.",
        ),
        IntegerField(
            "storage_used_bytes",
            required=False,
            disabled=True,
            label="Storage Used (Bytes)",
            help_text="Current storage usage for this project in bytes.",
        ),
        BooleanField("is_active", required=True, label="Is Active", help_text="Mark the project as active or not."),
    ]

    exclude_fields_from_list = ["description", "storage_used_bytes"]
    exclude_fields_from_create = ["storage_used_bytes", "is_active"]
    exclude_fields_from_edit = ["storage_used_bytes"]
    exclude_fields_from_detail = []

    def can_delete(self, request: Request) -> bool:
        """Disable deletion of projects. Projects can be soft deleted instead."""
        return False


class ProjectMembershipView(ModelView):
    """
    Custom admin panel View for the ProjectMembershipDB model.

    This view allows admins to manage project memberships: So assign users to projects
    and define their roles within the project.
    """

    page_size_options = PAGINATION_DEFAULTS
    exclude_fields_from_list = []
    exclude_fields_from_create = ["id", "created_at", "updated_at"]
    exclude_fields_from_edit = ["id", "created_at", "updated_at"]
    exclude_fields_from_detail = []


def register_admin_panel(app: FastAPI, engine: AsyncEngine) -> None:
    """
    Create and register an admin panel for the FastAPI app.
    """
    admin = Admin(engine, title="DivBase Admin")

    admin.add_view(UserView(UserDB, icon="fas fa-user", label="Users"))
    admin.add_view(ProjectView(ProjectDB, icon="fas fa-folder", label="Projects"))
    admin.add_view(ProjectMembershipView(ProjectMembershipDB, icon="fas fa-link", label="Project Memberships"))

    admin.mount_to(app)
