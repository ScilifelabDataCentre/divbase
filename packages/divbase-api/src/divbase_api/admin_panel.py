"""
Manage the starlette-admin interface for the DivBase API.

The views created for each model rely on overriding some of the default behavior provided by starlette-admin.
These overrides are on methods inside BaseModelView (parent of ModelView, which each custom class is inheriting from).
"""

import json
import logging
import pickle
from datetime import datetime, timezone
from typing import Any
from zoneinfo import ZoneInfo

from fastapi import FastAPI, Response
from pydantic import SecretStr
from sqlalchemy.exc import IntegrityError
from sqlalchemy.ext.asyncio import AsyncEngine
from starlette.requests import Request
from starlette_admin import (
    BaseAdmin,
    BooleanField,
    DateTimeField,
    EmailField,
    EnumField,
    HasOne,
    IntegerField,
    JSONField,
    StringField,
    TextAreaField,
)
from starlette_admin._types import RequestAction
from starlette_admin.auth import AdminUser, AuthProvider
from starlette_admin.contrib.sqla import Admin, ModelView
from starlette_admin.exceptions import FormValidationError

from divbase_api.db import get_db
from divbase_api.deps import _authenticate_frontend_user_from_tokens
from divbase_api.frontend_routes.auth import get_login, post_logout
from divbase_api.models.annoucements import AnnouncementDB, AnnouncementLevel, AnnouncementTarget
from divbase_api.models.project_versions import ProjectVersionDB
from divbase_api.models.projects import ProjectDB, ProjectMembershipDB
from divbase_api.models.revoked_tokens import RevokedTokenDB
from divbase_api.models.task_history import CeleryTaskMeta, TaskHistoryDB, TaskStartedAtDB
from divbase_api.models.users import UserDB
from divbase_api.security import get_password_hash

logger = logging.getLogger(__name__)

PAGINATION_DEFAULTS = [5, 10, 25, -1]  # (for number of items per page toggle)


def _format_cet_datetime(value: Any, field: Any, field_names: list[str]) -> str | None:
    """
    Helper function that can be called by overiding 'async def serialize_field_value()'
    for views that display datetimes.
    To change timezone, use patterns like this:
    dt = value.astimezone(timezone.utc) - displays UTC
    OR
    from zoneinfo import ZoneInfo
    dt = value.astimezone(ZoneInfo("Europe/Stockholm")) - displays CET
    """
    if isinstance(value, datetime) and field.name in field_names:
        dt = value.astimezone(ZoneInfo("Europe/Stockholm"))
        return dt.strftime("%Y-%m-%d %H:%M:%S %Z")
    return None


class UserView(ModelView):
    """
    Custom admin panel View for the UserDB model.

    Some Notes:
    - To handle password hashing, we override the create method to hash the password before calling the original create method.
    - Deletion of users is disabled, you can still mark a user as soft deleted though.
    - Adding/editing project memberships is disabled in this view, they should be handled in the ProjectMembership view.
    - Changing password for existing user is not supported: Users should use the email password reset flow instead.
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
        BooleanField("email_verified", help_text="Has the user verified their email address?"),
        DateTimeField(
            "last_password_change",
            help_text="Timestamp when the user last changed their password.",
            disabled=True,
        ),
        DateTimeField(
            "date_deleted",
            help_text="Timestamp when the user was soft deleted (else None). Value determined by system, cannot be edited.",
            disabled=True,
        ),
        "project_memberships",
        DateTimeField(
            "created_at", help_text="Timestamp when the entry was created. Value determined by system.", disabled=True
        ),
        DateTimeField(
            "updated_at",
            help_text="Timestamp when the entry was last updated. Value determined by system.",
            disabled=True,
        ),
    ]

    exclude_fields_from_list = ["hashed_password", "password"]
    exclude_fields_from_create = [
        "project_memberships",
        "is_deleted",
        "is_active",
        "last_password_change",
        "date_deleted",
    ]  # hashed_password wont be defined yet but needs to be included.
    exclude_fields_from_edit = ["project_memberships", "password", "hashed_password", "last_password_change"]
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

        return await super().create(request=request, data=data)

    async def edit(self, request: Request, pk: Any, data: dict) -> Any:
        """
        Override the default edit method to ensure that the `date_deleted` field is updated
        when/if a users soft deletion status is changed.
        """
        if "is_deleted" in data:
            if data["is_deleted"]:
                data["date_deleted"] = datetime.now(tz=timezone.utc)
            else:
                data["date_deleted"] = None

        return await super().edit(request=request, pk=pk, data=data)

    async def validate(self, request: Request, data: dict[str, Any]) -> None:
        """Custom validation to ensure a user cannot be both active and deleted at the same time."""

        # creation route does not set is_active/is_deleted fields, so we skip this check there
        # (otherwise keyerror)
        if "is_active" not in data or "is_deleted" not in data:
            return await super().validate(request=request, data=data)

        if data["is_active"] and data["is_deleted"]:
            raise FormValidationError(errors={"is_active": "Cannot set a user to active that is also set as deleted."})
        return await super().validate(request=request, data=data)

    def handle_exception(self, exc: Exception) -> None:
        """
        If an admin tries to create a user with an email that already exists, sqlalchemy will raise an IntegrityError.

        Handles catching and displaying this to avoid 500 internal server errors.
        """
        if isinstance(exc, IntegrityError):
            raise FormValidationError(errors={"email": "A user with this email already exists"})
        return super().handle_exception(exc)

    async def serialize_field_value(self, value: Any, field: Any, action: RequestAction, request: Request) -> Any:
        formatted = _format_cet_datetime(value, field, ["last_password_change", "date_deleted"])
        if formatted is not None:
            return formatted
        return await super().serialize_field_value(value, field, action, request)


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
        DateTimeField(
            "created_at", help_text="Timestamp when the entry was created. Value determined by system.", disabled=True
        ),
        DateTimeField(
            "updated_at",
            help_text="Timestamp when the entry was last updated. Value determined by system.",
            disabled=True,
        ),
    ]

    exclude_fields_from_list = ["description", "storage_used_bytes"]
    exclude_fields_from_create = ["storage_used_bytes", "is_active"]
    exclude_fields_from_edit = ["storage_used_bytes"]
    exclude_fields_from_detail = []

    def can_delete(self, request: Request) -> bool:
        """Disable deletion of projects. Projects can be soft deleted instead."""
        return False

    async def validate(self, request: Request, data: dict[str, Any]) -> None:
        """Custom validation to ensure a project name cannot have spaces and bucket name follows S3 bucket name rules."""

        if " " in data["name"]:
            raise FormValidationError(errors={"name": "Project name cannot contain spaces."})
        if " " in data["bucket_name"]:
            raise FormValidationError(errors={"bucket_name": "Bucket name cannot contain spaces."})

        if len(data["bucket_name"]) < 3 or len(data["bucket_name"]) > 63:
            raise FormValidationError(errors={"bucket_name": "Bucket name must be between 3 and 63 characters long."})

        return await super().validate(request=request, data=data)

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
    fields = [
        IntegerField("id", label="ID", disabled=True),
        HasOne("user", identity="user", label="User"),
        HasOne("project", identity="project", label="Project"),
        StringField("role", label="Role", required=True),
        DateTimeField("created_at", label="Created At", disabled=True),
        DateTimeField("updated_at", label="Updated At", disabled=True),
    ]

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

    async def serialize_field_value(self, value: Any, field: Any, action: RequestAction, request: Request) -> Any:
        formatted = _format_cet_datetime(value, field, ["created_at", "updated_at"])
        if formatted is not None:
            return formatted
        return await super().serialize_field_value(value, field, action, request)


class ProjectVersionsView(ModelView):
    """
    Custom admin panel View for the ProjectVersionDB model.
    """

    fields = [
        "id",
        StringField("name", required=True, label="Version Name", help_text="Unique name for the version."),
        TextAreaField("description", required=False, label="Description"),
        HasOne("project", identity="project", label="Project"),
        IntegerField(
            "user_id", label="User ID"
        ),  # No relationship created for this field in db model as this is for auditing only (can be null if user deleted)
        BooleanField("is_deleted", required=True, label="Is Deleted", help_text="Mark the version as deleted or not."),
        DateTimeField(
            "date_deleted",
            help_text="Timestamp when the version was soft deleted (else None). Value determined by system, cannot be edited.",
            disabled=True,
        ),
        JSONField(
            "files", required=True, label="Files", help_text="Mapping of file names to version IDs.", disabled=True
        ),
        DateTimeField(
            "created_at", help_text="Timestamp when the entry was created. Value determined by system.", disabled=True
        ),
        DateTimeField(
            "updated_at",
            help_text="Timestamp when the entry was last updated. Value determined by system.",
            disabled=True,
        ),
    ]

    page_size_options = PAGINATION_DEFAULTS
    exclude_fields_from_list = ["files"]
    exclude_fields_from_edit = ["id", "created_at", "updated_at", "files", "project", "user_id"]
    exclude_fields_from_detail = []

    def can_delete(self, request: Request) -> bool:
        """Disable deletion of project versions. Project versions can be soft deleted instead."""
        return False

    def can_create(self, request: Request) -> bool:
        """Disable creation of project versions. This is something users can create instead."""
        return False

    async def edit(self, request: Request, pk: Any, data: dict) -> Any:
        """
        Override the default edit method to ensure that the `date_deleted` field is updated
        when/if a version's soft deletion status is changed.
        """
        if "is_deleted" in data:
            if data["is_deleted"]:
                data["date_deleted"] = datetime.now(tz=timezone.utc)
            else:
                data["date_deleted"] = None

        return await super().edit(request=request, pk=pk, data=data)


class RevokedTokenView(ModelView):
    """
    Custom admin panel View for the RevokedTokenDB model.
    """

    page_size_options = PAGINATION_DEFAULTS
    fields = [
        IntegerField("id", label="ID", disabled=True),
        StringField("token_jti", label="Token JTI", required=True),
        StringField("token_type", label="Token Type", required=True),
        DateTimeField("revoked_at", label="Revoked At", disabled=True),
        StringField("revoked_reason", label="Revoked Reason", required=True),
        IntegerField("user_id", label="User ID", required=False),
        HasOne("user", identity="user", label="User"),
        DateTimeField("created_at", label="Created At", disabled=True),
        DateTimeField("updated_at", label="Updated At", disabled=True),
    ]

    exclude_fields_from_list = ["user_id"]
    exclude_fields_from_create = ["id", "user_id", "created_at", "updated_at", "revoked_at"]
    exclude_fields_from_edit = ["id", "user_id", "created_at", "updated_at", "revoked_reason"]
    exclude_fields_from_detail = []

    def handle_exception(self, exc: Exception) -> None:
        """
        Handles gracefully attempts to create/edit a revoked token entry that would otherwise become a 500 error:
            - A token_type that is not allowed (only refresh and password reset tokens can be revoked).
            - violating unique constraint on token_jti
        """
        if isinstance(exc, ValueError):  # raised by db models validate_token_type method
            raise FormValidationError(errors={"token_type": str(exc)})

        if isinstance(exc, IntegrityError):
            raise FormValidationError(errors={"token_jti": "A revoked token with this token_jti already exists."})

        return super().handle_exception(exc)

    async def serialize_field_value(self, value: Any, field: Any, action: RequestAction, request: Request) -> Any:
        formatted = _format_cet_datetime(value, field, ["created_at", "updated_at"])
        if formatted is not None:
            return formatted
        return await super().serialize_field_value(value, field, action, request)


class TaskHistoryView(ModelView):
    page_size_options = PAGINATION_DEFAULTS

    fields = [
        IntegerField("id", label="ID", disabled=True),
        StringField("task_id"),
        HasOne("user", identity="user", label="User"),
        HasOne("project", identity="project", label="Project"),
        HasOne("celery_meta", identity="celery-meta", label="Celery Task Details"),
        DateTimeField("created_at"),
    ]

    fields_default_sort = [("id", True)]  # False = descending, True = ascending

    async def serialize_field_value(self, value: Any, field: Any, action: RequestAction, request: Request) -> Any:
        """
        Override to format how values are displayed in the view.
        """
        if isinstance(value, datetime) and field.name in ["created_at"]:
            formatted = _format_cet_datetime(value, field, ["created_at"])
            if formatted is not None:
                return formatted
        return await super().serialize_field_value(value, field, action, request)

    def can_create(self, request: Request) -> bool:
        """Disable manual creation of task history entries."""
        return False

    def can_edit(self, request: Request) -> bool:
        """Optionally disable editing if task history should be read-only."""
        return False

    def can_delete(self, request: Request) -> bool:
        """Optionally disable deletion if task history should be immutable."""
        return False


class CeleryTaskMetaView(ModelView):
    """
    Custom admin panel View for CeleryTaskMeta (Celery's results backend table).

    Needs to be defined in order to display the entries in the admin panel, but is
    intended to only be viewed as a child of TaskHistoryView
    """

    page_size_options = PAGINATION_DEFAULTS
    exclude_fields_from_list = ["args", "kwargs", "result", "traceback"]

    fields = [
        IntegerField("id", label="ID", disabled=True),
        StringField("task_id", label="Task UUID", disabled=True),
        StringField("status", label="Celery Status", disabled=True),
        StringField("name", label="Task Name", disabled=True),
        StringField("worker", label="Worker", disabled=True),
        StringField("queue", label="Queue", disabled=True),
        IntegerField("retries", label="Retries", disabled=True),
        DateTimeField("date_done", label="Date Done", disabled=True),
        TextAreaField("args", label="Args", disabled=True),
        TextAreaField("kwargs", label="Kwargs", disabled=True),
        TextAreaField("result", label="Result", disabled=True),
        TextAreaField("traceback", label="Traceback", disabled=True),
    ]
    fields_default_sort = [("id", True)]  # False = descending, True = ascending

    async def serialize_field_value(self, value: Any, field: Any, action: RequestAction, request: Request) -> Any:
        """
        Override to deserialize Celery's binary fields for display.

        NOTE! serialize_field_value is a function in starlette-admin, so for the override to work, it cannot be renamed
        """
        if isinstance(value, datetime) and field.name == "date_done":
            formatted = _format_cet_datetime(value, field, ["date_done"])
            if formatted is not None:
                return formatted

        # For non-bytes values or fields we don't need to deserialize, use default behavior
        if not isinstance(value, bytes) or field.name not in ["args", "kwargs", "result"]:
            return await super().serialize_field_value(value, field, action, request)

        # NOTE: This is somewhat duplicated logic (also found in '_deserialize_celery_task_metadata' function from the task_history service layer).
        # It cannot be reused here as pydantic will raise validation errors as this function works on a per "cell" basis.
        # So the pydantic model cannot be created as would be missing all other (required) attributes.
        try:
            if field.name in ["args", "kwargs"]:
                # Deserialize args as JSON
                str_decoded = value.decode("utf-8") if isinstance(value, bytes) else value
                json_data = json.loads(str_decoded)
                return json.dumps(json_data, indent=2, default=str)

            elif field.name == "result":
                # Handle pickled results
                if isinstance(value, bytes) and value[:1] == b"\x80":
                    result_data = pickle.loads(value)
                else:
                    result_str = value.decode("utf-8") if isinstance(value, bytes) else value
                    result_data = json.loads(result_str)

                return json.dumps(result_data, indent=2, default=str)

        except Exception as e:
            logger.warning(f"Failed to deserialize {field.name}: {e}")
            return f"<Binary data: {len(value)} bytes>"

    def can_create(self, request: Request) -> bool:
        """Disable manual creation."""
        return False

    def can_edit(self, request: Request) -> bool:
        """Disable editing."""
        return False

    def can_delete(self, request: Request) -> bool:
        """Disable deletion."""
        return False


class TaskStartedAtView(ModelView):
    """
    Custom admin panel View for TaskStartedAtDB.
    """

    page_size_options = PAGINATION_DEFAULTS
    exclude_fields_from_list = ["created_at", "updated_at"]

    fields = [
        IntegerField("id", label="ID", disabled=True),
        StringField("task_id"),
        DateTimeField("started_at"),
    ]

    fields_default_sort = [("id", True)]  # False = descending, True = ascending

    async def serialize_field_value(self, value: Any, field: Any, action: RequestAction, request: Request) -> Any:
        """
        Override to format how values are displayed in the view.
        """
        if isinstance(value, datetime) and field.name in ["started_at"]:
            formatted = _format_cet_datetime(value, field, ["started_at"])
            if formatted is not None:
                return formatted
        return await super().serialize_field_value(value, field, action, request)

    def can_create(self, request: Request) -> bool:
        """Disable manual creation of task history entries."""
        return False

    def can_edit(self, request: Request) -> bool:
        """Optionally disable editing if task history should be read-only."""
        return False

    def can_delete(self, request: Request) -> bool:
        """Optionally disable deletion if task history should be immutable."""
        return False


class AnnouncementView(ModelView):
    """
    Custom admin panel View for the AnnouncementDB model.
    """

    page_size_options = PAGINATION_DEFAULTS
    fields = [
        IntegerField("id", label="ID", disabled=True),
        StringField("heading", label="Heading", required=True),
        TextAreaField("message", label="Message", required=False),
        EnumField("target", label="Target", required=True, enum=AnnouncementTarget),
        EnumField("level", label="Level", required=True, enum=AnnouncementLevel),
        DateTimeField(
            "auto_expire_at",
            help_text="Optional timestamp when the announcement will auto expire. After this timestamp the announcement will no longer be displayed on the frontend or CLI. Value can be left empty.",
        ),
        DateTimeField(
            "created_at", help_text="Timestamp when the entry was created. Value determined by system.", disabled=True
        ),
        DateTimeField(
            "updated_at",
            help_text="Timestamp when the entry was last updated. Value determined by system.",
            disabled=True,
        ),
    ]

    exclude_fields_from_list = ["message"]
    exclude_fields_from_edit = ["id", "created_at", "updated_at"]
    exclude_fields_from_detail = []


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
        # can't rely on dependency injection here like in FastAPI, so we manually obtain a db session
        async for db in get_db():
            logout_response = await post_logout(request, db)
        return logout_response

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

    admin.add_view(UserView(UserDB, icon="fas fa-user", label="Users", identity="user"))
    admin.add_view(ProjectView(ProjectDB, icon="fas fa-folder", label="Projects", identity="project"))
    admin.add_view(ProjectMembershipView(ProjectMembershipDB, icon="fas fa-link", label="Project Memberships"))
    admin.add_view(ProjectVersionsView(ProjectVersionDB, icon="fas fa-history", label="Project Versions"))
    admin.add_view(RevokedTokenView(RevokedTokenDB, icon="fas fa-ban", label="Revoked Tokens"))
    admin.add_view(TaskHistoryView(TaskHistoryDB, icon="fas fa-history", label="Task History"))
    admin.add_view(
        CeleryTaskMetaView(CeleryTaskMeta, icon="fas fa-tasks", label="Celery Task Meta", identity="celery-meta")
    )
    admin.add_view(TaskStartedAtView(TaskStartedAtDB, icon="fas fa-clock", label="Task Started At"))
    admin.add_view(
        AnnouncementView(AnnouncementDB, icon="fas fa-bullhorn", label="Announcements", identity="announcement")
    )
    admin.mount_to(app)
