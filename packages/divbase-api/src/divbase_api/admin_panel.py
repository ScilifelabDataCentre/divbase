"""
Manage the starlette-admin interface for the DivBase API.

The views created for each model rely on overriding some of the default behavior provided by starlette-admin.
These overrides are on methods inside BaseModelView (parent of ModelView, which each custom class is inheriting from).
"""

import logging
from typing import Any

from fastapi import FastAPI, Response
from pydantic import SecretStr
from sqlalchemy.exc import IntegrityError
from sqlalchemy.ext.asyncio import AsyncEngine
from starlette.requests import Request
from starlette_admin import BaseAdmin, BooleanField, EmailField, IntegerField, StringField, TextAreaField
from starlette_admin.auth import AdminUser, AuthProvider
from starlette_admin.contrib.sqla import Admin, ModelView
from starlette_admin.exceptions import FormValidationError

from divbase_api.db import get_db
from divbase_api.deps import _authenticate_frontend_user_from_tokens
from divbase_api.frontend_routes.auth import get_login, post_logout
from divbase_api.models.projects import ProjectDB, ProjectMembershipDB
from divbase_api.models.users import UserDB
from divbase_api.security import get_password_hash

logger = logging.getLogger(__name__)

PAGINATION_DEFAULTS = [5, 10, 25, -1]  # (for number of items per page toggle)


class UserView(ModelView):
    """
    Custom admin panel View for the UserDB model.

    Some Notes:
    - To handle password hashing, we override the create method to hash the password before calling the original create method.
    - Deletion of users is disabled, you can still mark a user as soft deleted though.
    - Adding/editing project memberships is disabled in this view, they should be handled in the ProjectMembership view.
    - Changing password for existing user is not supported: We will support a custom password reset flow instead.
    """

    page_size_options = PAGINATION_DEFAULTS
    fields = [
        "id",
        StringField("name", required=True, help_text="Full name of the user."),
        EmailField("email", required=True, help_text="Email address of the user."),
        StringField("password", required=True, help_text="Password for the user."),
        StringField(
            "hashed_password",
            required=False,
            disabled=True,
            help_text="Hashed password is auto created by the system from password, cannot be edited.",
        ),
        BooleanField("is_admin", help_text="Is the user an admin?"),
        BooleanField("is_active", help_text="Is the user active?"),
        BooleanField("is_deleted", help_text="Is the user deleted?"),
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
        and create a hashed password which will be needed to add to the db before calling the original create method.
        """
        password = data["password"]
        hashed_password = get_password_hash(SecretStr(password))
        data["hashed_password"] = hashed_password

        return await super().create(request, data)

    async def validate(self, request: Request, data: dict[str, Any]) -> None:
        """Custom validation to ensure a user cannot be both active and deleted at the same time."""

        # creation route does not set is_active/is_deleted fields, so we skip this check there
        # (otherwise keyerror)
        if "is_active" not in data or "is_deleted" not in data:
            return await super().validate(request, data)

        if data["is_active"] and data["is_deleted"]:
            raise FormValidationError(errors={"is_active": "Cannot set a user to active that is also set as deleted."})
        return await super().validate(request, data)

    def handle_exception(self, exc: Exception) -> None:
        """
        If an admin tries to create a user with an email that already exists, sqlalchemy will raise an IntegrityError.

        Handles catching and displaying this to avoid 500 internal server errors.
        """
        if isinstance(exc, IntegrityError):
            raise FormValidationError(errors={"email": "A user with this email already exists"})
        return super().handle_exception(exc)


class ProjectView(ModelView):
    """
    Custom admin panel View for the ProjectDB model.

    Project memberships are managed in the ProjectMembership view.
    """

    page_size_options = PAGINATION_DEFAULTS
    fields = [
        "id",
        StringField(
            "name", required=True, label="Project Name", help_text="Unique name for the project, no spaces allowed."
        ),
        TextAreaField("description", required=False, label="Description"),
        StringField(
            "bucket_name", required=True, label="Bucket Name", help_text="Unique S3 bucket name for the project."
        ),
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

    def handle_exception(self, exc: Exception) -> None:
        """
        sqlalchemy will raise an IntegrityError if an admin tries to create/edit a project to have either a
         (1) project name that already exists or
         (2) bucket name that already exists.

        Handles catching and displaying this to avoid 500 internal server errors.
        """
        if isinstance(exc, IntegrityError):
            orig_message = str(exc.orig).lower() if exc.orig else str(exc).lower()

            # As multiple fields could have raised the integrity error, we need to see which one it was.
            if "bucket_name" in orig_message:
                raise FormValidationError(errors={"bucket_name": "A project with this bucket name already exists."})
            if "name" in orig_message:
                raise FormValidationError(errors={"name": "A project with this name already exists."})

        return super().handle_exception(exc)


class ProjectMembershipView(ModelView):
    """
    Custom admin panel View for the ProjectMembershipDB model.

    This view allows admins to manage project memberships: So assign users to projects
    and (re)define their roles within the project.
    """

    page_size_options = PAGINATION_DEFAULTS
    exclude_fields_from_list = []
    exclude_fields_from_create = ["id", "created_at", "updated_at"]
    exclude_fields_from_edit = ["id", "user_id", "project_id", "created_at", "updated_at"]
    exclude_fields_from_detail = []

    def handle_exception(self, exc: Exception) -> None:
        """
        sqlalchemy will raise an IntegrityError if an admin tries to create a new projectmembership if one already exists
        for the given user id + project id combo.

        Handles catching and displaying this to avoid a 500 internal server error.
        """
        if isinstance(exc, IntegrityError):
            raise FormValidationError(
                errors={
                    "role": """A project membership table entry with this user id + project id already exists. 
                    Edit that entry instead or delete the entry first."""
                }
            )

        return super().handle_exception(exc)


class DivBaseAuthProvider(AuthProvider):
    """
    This class enables starlette-admin to make use of DivBase's pre-existing auth system.

    The methods below are overriding several existing methods in the AuthProvider class (and its parent BaseAuthProvider).
    """

    async def render_login(self, request: Request, admin: BaseAdmin) -> Response:
        """Override the default starlette-admin login method to use our frontend get_login route/page."""
        return await get_login(request)

    async def render_logout(self, request: Request, admin: BaseAdmin) -> Response:
        """Override the default starlette-admin logout to use our frontend post_logout function/route."""
        return await post_logout(request)

    async def is_authenticated(self, request: Request) -> bool:
        """
        Overrides the default implementation to use our pre-existing DivBase auth system.

        In which we determine the user from their JWT tokens stored in httponly cookies.

        As expected, we also confirm user.is_admin for access.
        """
        access_token = request.cookies.get("access_token")
        refresh_token = request.cookies.get("refresh_token")
        if not access_token and not refresh_token:
            return False

        try:
            # Starlette does not support dependency injection like FastAPI,
            # so we need to manually obtain the database session here.
            async for db in get_db():
                user = await _authenticate_frontend_user_from_tokens(
                    access_token=access_token, refresh_token=refresh_token, db=db
                )

                if user and user.is_admin and user.is_active:
                    # Store user info in the request state so it can be accessed by e.g. get_admin_user
                    request.state.user = {"id": user.id, "name": user.name, "is_admin": user.is_admin}
                    return True
        except Exception as e:
            logger.warning(
                f"An error occurred while attempting to authenticate a user on the starlette-admin panel, details: {e}"
            )
            return False

        return False

    def get_admin_user(self, request: Request) -> AdminUser | None:
        """
        Retrieve the current (admin) user for display on the admin panel.

        This controls the display of the user info in the top right of the admin panel and makes having a logout button possible.
        """
        user = request.state.user
        if not user:
            return None

        return AdminUser(username=user["name"], photo_url=None)


def register_admin_panel(app: FastAPI, engine: AsyncEngine) -> None:
    """
    Create and register an admin panel for the FastAPI app.
    """
    admin = Admin(engine=engine, title="DivBase Admin", auth_provider=DivBaseAuthProvider())

    admin.add_view(UserView(UserDB, icon="fas fa-user", label="Users"))
    admin.add_view(ProjectView(ProjectDB, icon="fas fa-folder", label="Projects"))
    admin.add_view(ProjectMembershipView(ProjectMembershipDB, icon="fas fa-link", label="Project Memberships"))

    admin.mount_to(app)
