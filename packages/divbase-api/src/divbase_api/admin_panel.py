"""
Manage the starlette-admin interface for the DivBase API.
"""

from typing import Any

from fastapi import FastAPI
from pydantic import SecretStr
from sqlalchemy.ext.asyncio import AsyncEngine
from starlette.requests import Request
from starlette_admin import EmailField, StringField
from starlette_admin.contrib.sqla import Admin, ModelView

from divbase_api.models.projects import ProjectDB, ProjectMembershipDB
from divbase_api.models.users import UserDB
from divbase_api.security import get_password_hash


class UserView(ModelView):
    """
    Custom admin panel View for the UserDB model.

    These settings override some of the defaults set by starlette_admin BaseModelView (parent of ModelView).

    Some Notes:
    - To handle password hashing, we override the create method to hash the password before calling the original create method.
    - We disable deletion of users, you can mark a user as soft deleted though.
    - Adding/editing project memberships for a user is disabled, they should be handled in the ProjectMembership view.
    - Changing password for existing user is not supported: TODO see if can fix or explain here.
    """

    page_size_options = [5, 10, 25, -1]  # (for number of items per page toggle)
    fields = [
        "id",
        "name",
        EmailField("email"),
        StringField("password", label="password", required=True),
        StringField(
            "hashed_password",
            label="hashed password",
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
    exclude_fields_from_create = ["project_memberships", "is_deleted", "is_active"]
    exclude_fields_from_edit = ["password", "hashed_password", "project_memberships"]
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


def register_admin_panel(app: FastAPI, engine: AsyncEngine) -> None:
    """
    Create and register an admin panel for the FastAPI app.
    """
    admin = Admin(engine, title="DivBase Admin")

    admin.add_view(UserView(UserDB, icon="fas fa-user"))
    admin.add_view(ModelView(ProjectDB, icon="fas fa-folder"))
    admin.add_view(ModelView(ProjectMembershipDB, icon="fas fa-link"))

    admin.mount_to(app)
