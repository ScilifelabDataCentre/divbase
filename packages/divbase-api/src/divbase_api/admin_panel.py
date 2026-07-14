"""
Manages the starlette-admin interface for the DivBase API.

### Inheritance structure of the views:
- The views created for each model rely on overriding some of the default behavior provided by starlette-admin.
- We have a DivBaseModelView class that all concrete models inherit from with some default behaviours and settings.
- DivBaseModelView in turn inherits from (starlette-admin's) ModelView which in turn inherits from BaseModelView.

### Datetime/timezone handling:
- We use the TimezoneConfig from starlette-admin to control datetime fields in the admin panel.
- All datetime fields are displayed in the admin panel in Europe/Stockholm timezone (CET/CEST).
- All datetime fields are stored in the database in UTC timezone (starlette-admin handles the conversion for us).
- To make clear to admins that the datetimes they are looking at are in Europe/Stockholm timezone,
we add a label to all datetime fields in the admin panel to indicate this.
"""

import json
import pickle
from datetime import datetime, timezone
from typing import Any

import structlog
from fastapi import FastAPI, Response
from pydantic import SecretStr
from sqlalchemy.exc import IntegrityError
from sqlalchemy.ext.asyncio import AsyncEngine
from starlette.requests import Request
from starlette_admin import (
    BaseAdmin,
    BaseField,
    BooleanField,
    DateTimeField,
    EmailField,
    EnumField,
    HasOne,
    IntegerField,
    JSONField,
    StringField,
    TextAreaField,
    TimezoneConfig,
)
from starlette_admin._types import RequestAction
from starlette_admin.auth import AdminUser, AuthProvider
from starlette_admin.contrib.sqla import Admin, ModelView
from starlette_admin.exceptions import FormValidationError

from divbase_api.db import get_db
from divbase_api.deps import _authenticate_frontend_user_from_tokens
from divbase_api.frontend_routes.auth import get_login, post_logout
from divbase_api.models.announcements import AnnouncementDB, AnnouncementLevel, AnnouncementTarget
from divbase_api.models.personal_access_tokens import PersonalAccessTokenDB
from divbase_api.models.project_versions import ProjectVersionDB
from divbase_api.models.projects import ProjectDB, ProjectMembershipDB, ProjectRoles
from divbase_api.models.queue_status import QueueStatusDB
from divbase_api.models.revoked_tokens import RevokedTokenDB, TokenRevokeReason
from divbase_api.models.task_history import CeleryTaskMeta, TaskHistoryDB, TaskStartedAtDB
from divbase_api.models.users import UserDB
from divbase_api.security import TokenType, get_password_hash
from divbase_lib.divbase_constants import DIVBASE_SERVER_TIMEZONE

logger = structlog.get_logger(__name__)

# We add this label to all DateTimeField instances to make it clear to the admin what timezone they are working in.
DATETIME_TIMEZONE_LABEL = f"({DIVBASE_SERVER_TIMEZONE} time)"


def _basedb_model_fields() -> list[BaseField]:
    """Helper fn to append the common BaseDBModel fields to a ModelView's field list."""
    id = IntegerField("id", label="ID", disabled=True)
    created_at = DateTimeField(
        "created_at",
        label=f"Created At {DATETIME_TIMEZONE_LABEL}",
        help_text="Timestamp when the entry was created. Value determined by system.",
        disabled=True,
    )
    updated_at = DateTimeField(
        "updated_at",
        label=f"Updated At {DATETIME_TIMEZONE_LABEL}",
        help_text="Timestamp when the entry was last updated. Value determined by system.",
        disabled=True,
    )
    return [id, created_at, updated_at]


def _is_deleted_date_deleted_fields() -> list[BaseField]:
    """Helper fn to append 'is_deleted' and 'date_deleted' fields to a ModelView's field list."""
    is_deleted = BooleanField(
        "is_deleted", required=True, label="Is Deleted", help_text="Mark the entry as soft deleted or not."
    )
    date_deleted = DateTimeField(
        "date_deleted",
        label=f"Date Deleted {DATETIME_TIMEZONE_LABEL}",
        help_text="Timestamp when the entry was soft deleted (else None). Value determined by system, cannot be edited.",
        disabled=True,
    )
    return [is_deleted, date_deleted]


class DivBaseModelView(ModelView):
    """Shared admin view for all DivBase DB models."""

    page_size_options = [5, 10, 25, -1]  # (for number of items per page toggle)
    fields_default_sort = [("id", True)]  # False = descending, True = ascending

    exclude_fields_from_list = ["created_at", "updated_at"]
    exclude_fields_from_create = ["id", "created_at", "updated_at"]
    exclude_fields_from_edit = ["id", "created_at", "updated_at"]
    exclude_fields_from_detail = []

    async def edit(self, request: Request, pk: Any, data: dict) -> Any:
        """
        Override the edit method to keep is_deleted and date_deleted in sync when editing is_deleted.
        Models without "is_deleted"/"date_deleted" fields are not impacted.
        """
        if "is_deleted" in data:
            if data["is_deleted"]:
                data["date_deleted"] = datetime.now(tz=timezone.utc)
            else:
                data["date_deleted"] = None
        return await super().edit(request=request, pk=pk, data=data)


class UserView(DivBaseModelView):
    """
    Custom admin panel View for the UserDB model.

    Some Notes:
    - To handle password hashing, we override the create method to hash the password before calling the original create method.
    - Deletion of users is disabled, you can still mark a user as soft deleted though.
    - Adding/editing project memberships is disabled in this view, they should be handled in the ProjectMembership view.
    - Changing password for existing user is not supported: Users should use the email password reset flow instead.
    """

    fields = (
        [
            StringField("name", required=True, help_text="Full name of the user."),
            EmailField("email", required=True, help_text="Email address of the user."),
            StringField("organisation", required=True, help_text="Organisation of the user."),
            StringField("organisation_role", required=True, help_text="Role of the user within their organisation."),
            StringField("password", required=True, help_text="Password for the user."),
            StringField(
                "hashed_password",
                required=False,
                disabled=True,
                help_text="Hashed password is auto created by the system from password, cannot be edited.",
            ),
            BooleanField("is_admin", help_text="Is the user an admin?"),
            BooleanField("is_active", help_text="Is the user active?"),
            BooleanField("email_verified", help_text="Has the user verified their email address?"),
            DateTimeField(
                "last_password_change",
                label=f"Last Password Change {DATETIME_TIMEZONE_LABEL}",
                help_text="Timestamp when the user last changed their password.",
                disabled=True,
            ),
            "project_memberships",
        ]
        + _is_deleted_date_deleted_fields()
        + _basedb_model_fields()
    )

    exclude_fields_from_list = [
        "hashed_password",
        "password",
        "organisation_role",
        "date_deleted",
        "last_password_change",
        "created_at",
        "updated_at",
    ]
    exclude_fields_from_create = [
        "id",
        "created_at",
        "updated_at",
        "project_memberships",
        "is_deleted",
        "is_active",
        "last_password_change",
        "date_deleted",
    ]
    exclude_fields_from_edit = [
        "id",
        "created_at",
        "updated_at",
        "project_memberships",
        "password",
        "hashed_password",  # hashed_password wont be defined yet but needs to be included here.
        "last_password_change",
    ]
    exclude_fields_from_detail = ["hashed_password", "password"]

    def can_delete(self, request: Request) -> bool:
        """Can only soft delete users"""
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

    async def validate(self, request: Request, data: dict[str, Any]) -> None:
        """Custom validation to ensure a user cannot be both active and deleted at the same time."""

        # creation route does not set is_active/is_deleted fields, so we skip this check there
        # (otherwise keyerror)
        if "is_active" not in data or "is_deleted" not in data:
            return await super().validate(request=request, data=data)

        if data["is_active"] and data["is_deleted"]:
            raise FormValidationError(errors={"is_active": "Cannot set a user as both active and deleted."})
        return await super().validate(request=request, data=data)

    def handle_exception(self, exc: Exception) -> None:
        """
        If an admin tries to create a user with an email that already exists, sqlalchemy will raise an IntegrityError.

        Handles catching and displaying this to avoid 500 internal server errors.
        """
        if isinstance(exc, IntegrityError):
            raise FormValidationError(errors={"email": "A user with this email already exists"})
        return super().handle_exception(exc)


class ProjectView(DivBaseModelView):
    """
    Custom admin panel View for the ProjectDB model.

    Project memberships are managed in the ProjectMembership view.
    """

    fields = (
        [
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
        + _is_deleted_date_deleted_fields()
        + _basedb_model_fields()
    )

    exclude_fields_from_list = ["created_at", "updated_at", "description", "storage_used_bytes"]
    exclude_fields_from_create = ["id", "created_at", "updated_at", "storage_used_bytes", "is_active", "is_deleted"]
    exclude_fields_from_edit = ["id", "created_at", "updated_at", "storage_used_bytes"]

    def can_delete(self, request: Request) -> bool:
        return False

    async def validate(self, request: Request, data: dict[str, Any]) -> None:
        """Custom validation to ensure a project name cannot have spaces and bucket name follows S3 bucket name rules."""

        if " " in data["name"]:
            raise FormValidationError(errors={"name": "Project name cannot contain spaces."})
        if " " in data["bucket_name"]:
            raise FormValidationError(errors={"bucket_name": "Bucket name cannot contain spaces."})

        if len(data["bucket_name"]) < 3 or len(data["bucket_name"]) > 63:
            raise FormValidationError(errors={"bucket_name": "Bucket name must be between 3 and 63 characters long."})

        # creation route does not set is_active/is_deleted fields, so we skip this check there
        # (otherwise keyerror)
        if "is_active" not in data or "is_deleted" not in data:
            return await super().validate(request=request, data=data)

        if data["is_active"] and data["is_deleted"]:
            raise FormValidationError(errors={"is_active": "Cannot set a project as both active and deleted."})
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


class ProjectMembershipView(DivBaseModelView):
    """
    Custom admin panel View for the ProjectMembershipDB model.

    This view allows admins to manage project memberships: So assign users to projects
    and (re)define their roles within the project.
    """

    fields = [
        HasOne("user", identity="user", label="User"),
        HasOne("project", identity="project", label="Project"),
        EnumField("role", label="Role", required=True, enum=ProjectRoles),
    ] + _basedb_model_fields()

    exclude_fields_from_list = []
    exclude_fields_from_edit = ["id", "created_at", "updated_at", "user_id", "project_id"]

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


class ProjectVersionsView(DivBaseModelView):
    """Custom admin panel View for the ProjectVersionDB model."""

    fields = (
        [
            StringField("name", required=True, label="Version Name", help_text="Unique name for the version."),
            TextAreaField("description", required=False, label="Description"),
            HasOne("project", identity="project", label="Project"),
            IntegerField(
                "user_id", label="User ID"
            ),  # No relationship created for this field in db model as this is for auditing only (can be null if user deleted)
            JSONField(
                "files", required=True, label="Files", help_text="Mapping of file names to version IDs.", disabled=True
            ),
        ]
        + _is_deleted_date_deleted_fields()
        + _basedb_model_fields()
    )

    exclude_fields_from_list = ["files"]
    exclude_fields_from_edit = ["id", "created_at", "updated_at", "files", "project", "user_id"]

    def can_delete(self, request: Request) -> bool:
        """Disable deletion of project versions. Project versions can be soft deleted instead."""
        return False

    def can_create(self, request: Request) -> bool:
        """Disable creation of project versions. This is something users can create instead."""
        return False


class RevokedTokenView(DivBaseModelView):
    """Custom admin panel View for the RevokedTokenDB model."""

    fields = [
        StringField("token_jti", label="Token JTI", required=True),
        EnumField("token_type", label="Token Type", required=True, enum=TokenType),
        DateTimeField("revoked_at", label=f"Revoked At {DATETIME_TIMEZONE_LABEL}", disabled=True),
        EnumField("revoked_reason", label="Revoke Reason", required=True, enum=TokenRevokeReason),
        IntegerField("user_id", label="User ID", required=False),
        HasOne("user", identity="user", label="User"),
    ] + _basedb_model_fields()

    exclude_fields_from_list = ["user_id"]
    exclude_fields_from_create = ["id", "created_at", "updated_at", "user_id", "revoked_at"]
    exclude_fields_from_edit = ["id", "created_at", "updated_at", "user_id", "revoked_at", "revoked_reason"]

    def handle_exception(self, exc: Exception) -> None:
        """
        Handles gracefully attempts to create/edit a revoked token entry that would otherwise become a 500 error:
            - A token_type that is not allowed (only refresh and password reset tokens can be revoked).
            - violating unique constraint on token_jti
        """
        if isinstance(exc, ValueError):  # raised by db models validate_token_type method
            raise FormValidationError(errors={"token_type": str(exc)})

        if isinstance(exc, IntegrityError):
            orig_message = str(exc.orig).lower() if exc.orig else str(exc).lower()

            if "token_jti" in orig_message and "unique constraint" in orig_message:
                raise FormValidationError(errors={"token_jti": "A revoked token with this token_jti already exists."})
            else:
                raise FormValidationError(
                    errors={"token_jti": f"Unexpected integrity error: {orig_message} " + orig_message}
                )

        return super().handle_exception(exc)


class TaskHistoryView(DivBaseModelView):
    fields = [
        StringField("task_id"),
        HasOne("user", identity="user", label="User"),
        HasOne("project", identity="project", label="Project"),
        HasOne("celery_meta", identity="celery-meta", label="Celery Task Details"),
    ] + _basedb_model_fields()

    exclude_fields_from_list = ["id", "updated_at"]

    def can_create(self, request: Request) -> bool:
        """task history entries should be immutable"""
        return False

    def can_edit(self, request: Request) -> bool:
        """task history entries should be immutable"""
        return False

    def can_delete(self, request: Request) -> bool:
        """task history entries should be immutable"""
        return False


class CeleryTaskMetaView(DivBaseModelView):
    """
    Custom admin panel View for CeleryTaskMeta (Celery's results backend table).

    Needs to be defined in order to display the entries in the admin panel, but is
    intended to only be viewed as a child of TaskHistoryView
    """

    fields = [
        IntegerField("id", label="ID", disabled=True),
        StringField("task_id", label="Task UUID", disabled=True),
        StringField("status", label="Celery Status", disabled=True),
        StringField("name", label="Task Name", disabled=True),
        StringField("worker", label="Worker", disabled=True),
        StringField("queue", label="Queue", disabled=True),
        IntegerField("retries", label="Retries", disabled=True),
        DateTimeField("date_done", label=f"Date Done {DATETIME_TIMEZONE_LABEL}", disabled=True),
        TextAreaField("args", label="Args", disabled=True),
        TextAreaField("kwargs", label="Kwargs", disabled=True),
        TextAreaField("result", label="Result", disabled=True),
        TextAreaField("traceback", label="Traceback", disabled=True),
    ]
    exclude_fields_from_list = ["args", "kwargs", "result", "traceback"]

    async def serialize_field_value(self, value: Any, field: Any, action: RequestAction, request: Request) -> Any:
        """Override to deserialize Celery's binary fields for display."""
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
        """Celery managed table, no edits or creation allowed."""
        return False

    def can_edit(self, request: Request) -> bool:
        """Celery managed table, no edits or creation allowed."""
        return False

    def can_delete(self, request: Request) -> bool:
        """Celery managed table, no edits or creation allowed."""
        return False


class TaskStartedAtView(DivBaseModelView):
    """Custom admin panel View for TaskStartedAtDB."""

    fields = [
        StringField("task_id"),
        DateTimeField("started_at", label=f"Started At {DATETIME_TIMEZONE_LABEL}", disabled=True),
    ] + _basedb_model_fields()

    def can_create(self, request: Request) -> bool:
        """System managed table, no edits or creation allowed."""
        return False

    def can_edit(self, request: Request) -> bool:
        """System managed table, no edits or creation allowed."""
        return False

    def can_delete(self, request: Request) -> bool:
        """System managed table, no edits or creation allowed."""
        return False


class AnnouncementView(DivBaseModelView):
    """Custom admin panel View for the AnnouncementDB model."""

    fields = [
        StringField("heading", label="Heading", required=True),
        TextAreaField("message", label="Message", required=False),
        EnumField("target", label="Target", required=True, enum=AnnouncementTarget),
        EnumField("level", label="Level", required=True, enum=AnnouncementLevel),
        DateTimeField(
            "auto_expire_at",
            label=f"Auto Expire At {DATETIME_TIMEZONE_LABEL}",
            help_text="Optional timestamp when the announcement will auto expire. After this timestamp the announcement will no longer be displayed on the frontend or CLI. Value can be left empty.",
        ),
    ] + _basedb_model_fields()

    exclude_fields_from_list = ["message"]


class QueueStatusView(DivBaseModelView):
    """Custom admin panel View for the QueueStatusDB model."""

    fields = [
        BooleanField(
            "is_closed",
            help_text=(
                "Is the queue closed for new tasks? "
                "If the scheduled start time is in the past or not set, the queue is closed. "
                "If scheduled start time is in the future, queue will be closed at that time."
            ),
        ),
        DateTimeField(
            "scheduled_start",
            label=f"Scheduled closure start time {DATETIME_TIMEZONE_LABEL}",
            help_text=(
                "Optional: When should the queue closure set above take effect? "
                "Leave empty for queue to be closed immediatialy "
                "If the queue is open and you set this field, it will do nothing."
            ),
        ),
        TextAreaField(
            "reason_for_users",
            help_text=(
                "Message shown to users when the queue is closed (max 500 characters). "
                "You could write: "
                "The queuing system is currently closed for new tasks due to planned upcoming maintenance. Please try again later."
            ),
        ),
    ] + _basedb_model_fields()

    def can_delete(self, request: Request) -> bool:
        """Disable deletion - this is a singleton table with only 1 row."""
        return False

    def can_create(self, request: Request) -> bool:
        """Disable creation - this is a singleton table with only 1 row (the row is created by the alembic migration)"""
        return False

    async def edit(self, request: Request, pk: Any, data: dict) -> Any:
        """keep `scheduled_start` field in sync with `is_closed` is set to False."""
        if not data.get("is_closed"):
            data["scheduled_start"] = None
        return await super().edit(request=request, pk=pk, data=data)

    async def validate(self, request: Request, data: dict[str, Any]) -> None:
        errors: dict[str, str] = {}

        if data.get("scheduled_start") and not data.get("is_closed"):
            errors["scheduled_start"] = "Cannot have a scheduled start time without the queue being closed."

        reason = data.get("reason_for_users") or ""
        if data.get("is_closed") and not reason.strip():
            errors["reason_for_users"] = "Reason for users is required when the queue is closed."
        if len(reason) > 500:
            errors["reason_for_users"] = "Message must be 500 characters or less."

        if errors:
            raise FormValidationError(errors=errors)
        return await super().validate(request=request, data=data)


class PersonalAccessTokenView(DivBaseModelView):
    """
    Custom admin panel View for the PersonalAccessTokenDB model.
    As with passwords, the hashed_token is never displayed in the admin panel on purpose.
    """

    fields = (
        [
            StringField("name", label="Name", required=True),
            TextAreaField("description", required=False, label="Description"),
            JSONField(
                "permissions",
                required=False,
                label="Permissions",
                help_text="Permissions associated with the PAT.",
                disabled=True,
            ),
            DateTimeField(
                "expires_at",
                label=f"Expires At {DATETIME_TIMEZONE_LABEL}",
                help_text="Timestamp when the PAT expires. Value can be left empty for no expiration. Value determined by system, cannot be edited.",
                required=False,
                disabled=True,
            ),
            DateTimeField(
                "last_used_at",
                label=f"Last Used At {DATETIME_TIMEZONE_LABEL}",
                help_text="Timestamp when the PAT was last used. Value empty if not yet used. Value determined by system, cannot be edited.",
                required=False,
                disabled=True,
            ),
            IntegerField("user_id", label="User ID", required=False),
            HasOne("user", identity="user", label="User"),
        ]
        + _is_deleted_date_deleted_fields()
        + _basedb_model_fields()
    )

    exclude_fields_from_list = ["hashed_token", "permissions", "user_id", "id", "created_at", "updated_at"]
    exclude_fields_from_edit = ["hashed_token", "id", "created_at", "updated_at", "user_id", "user", "permissions"]
    exclude_fields_from_detail = ["hashed_token", "user_id"]

    def can_create(self, request: Request) -> bool:
        """Should create these in the frontend UI instead"""
        return False


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
    """Create and register an admin panel for the FastAPI app."""
    timezone_config = TimezoneConfig(
        default_timezone=DIVBASE_SERVER_TIMEZONE,
        database_timezone="UTC",
        timezone_cookie_name=None,
    )
    admin = Admin(
        engine=engine,
        title="DivBase Admin",
        auth_provider=DivBaseAuthProvider(),
        timezone_config=timezone_config,
    )

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
    admin.add_view(
        QueueStatusView(QueueStatusDB, icon="fas fa-power-off", label="Queue Status", identity="queue-status")
    )
    admin.add_view(
        PersonalAccessTokenView(
            PersonalAccessTokenDB, icon="fas fa-key", label="Personal Access Tokens", identity="personal-access-token"
        )
    )
    admin.mount_to(app)
