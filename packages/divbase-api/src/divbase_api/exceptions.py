"""
DivBase API custom exceptions.

All exceptions in this module should inherit from DivBaseAPIException.
"""

import logging
from pathlib import Path

from fastapi import status


class DivBaseAPIException(Exception):
    """
    Base exception for all DivBaseAPI errors

    Each subclass declares how `exception_handlers.divbase_api_exception_handler` should treat it:
    - log_level: What log level to use
    - include_exc_info: whether to include the exception traceback in the logs.
    - frontend_redirect_url: if set, frontend (non-API) requests are redirected here instead of shown the
      generic error page.
    - frontend_message: if set, frontend requests see this instead of `message` (API responses always see
      `message`) - used when the raw message is too internal/technical for a human-facing page.
    """

    log_level: int = logging.INFO
    include_exc_info: bool = False
    frontend_redirect_url: str | None = None
    frontend_message: str | None = None

    def __init__(self, message: str, status_code: int, headers: dict[str, str] | None = None):
        self.message = message
        self.status_code = status_code
        self.headers = headers or {}
        super().__init__(self.message)

    @property
    def error_type(self) -> str:
        """Returned in API responses as the (error) type."""
        return type(self).__name__


class AuthenticationError(DivBaseAPIException):
    log_level = logging.INFO
    include_exc_info = False
    frontend_redirect_url = "/login"

    def __init__(self, message: str = "Authentication required"):
        default_headers = {"WWW-Authenticate": "Bearer"}
        super().__init__(message=message, status_code=status.HTTP_401_UNAUTHORIZED, headers=default_headers)


class AuthorizationError(DivBaseAPIException):
    log_level = logging.INFO
    include_exc_info = False
    frontend_redirect_url = "/login"

    def __init__(self, message: str = "Authorization required"):
        default_headers = {"WWW-Authenticate": "Bearer"}
        super().__init__(message=message, status_code=status.HTTP_403_FORBIDDEN, headers=default_headers)


class UserRegistrationError(DivBaseAPIException):
    """Exception for failed registration of new user or update of existing user."""

    log_level = logging.INFO
    include_exc_info = False

    def __init__(
        self,
        internal_logging_message: str,  # for logs, can be more detailed
        user_message: str = "Registration failed. Please try again.",  # given to end user, generic
    ):
        self.user_message = user_message
        super().__init__(message=internal_logging_message, status_code=status.HTTP_400_BAD_REQUEST)


class ProjectNotFoundError(DivBaseAPIException):
    log_level = logging.DEBUG
    include_exc_info = False
    frontend_message = "Project not found or you don't have access."

    def __init__(self, message: str = "Project not found or you don't have access"):
        super().__init__(message=message, status_code=status.HTTP_404_NOT_FOUND)


class ProjectMemberNotFoundError(DivBaseAPIException):
    log_level = logging.INFO
    include_exc_info = False
    frontend_redirect_url = "/projects"

    def __init__(self, message: str = "Project member not found"):
        super().__init__(message=message, status_code=status.HTTP_404_NOT_FOUND)


class UserNotFoundError(DivBaseAPIException):
    log_level = logging.DEBUG
    include_exc_info = False

    def __init__(self, message: str = "User not found"):
        super().__init__(message=message, status_code=status.HTTP_404_NOT_FOUND)


class ProjectMemberAlreadyExistsError(DivBaseAPIException):
    log_level = logging.DEBUG
    include_exc_info = False

    def __init__(self, message: str = "User is already a member of this project"):
        super().__init__(message=message, status_code=status.HTTP_409_CONFLICT)


class ProjectCreationError(DivBaseAPIException):
    log_level = logging.WARNING
    include_exc_info = True

    def __init__(self, message: str = "Project creation failed"):
        super().__init__(message=message, status_code=status.HTTP_400_BAD_REQUEST)


class ProjectVersionCreationError(DivBaseAPIException):
    """Raised when there is an error creating a new project version."""

    log_level = logging.INFO
    include_exc_info = False

    def __init__(self, message: str = "Failed to create a new project version"):
        super().__init__(message=message, status_code=status.HTTP_400_BAD_REQUEST)


class ProjectVersionAlreadyExistsError(DivBaseAPIException):
    """Raised when attempting to create a project version that already exists."""

    log_level = logging.INFO
    include_exc_info = False

    def __init__(
        self, message: str = "A project version with this name already exists, please choose a different name."
    ):
        super().__init__(message=message, status_code=status.HTTP_400_BAD_REQUEST)


class ProjectVersionNotFoundError(DivBaseAPIException):
    """Raised when a project version is not found."""

    log_level = logging.INFO
    include_exc_info = False

    def __init__(self, message: str = "Could not find the specified project version"):
        super().__init__(message=message, status_code=status.HTTP_404_NOT_FOUND)


class ProjectVersionSoftDeletedError(DivBaseAPIException):
    """Raised when a user tries to modify a project version that is soft-deleted."""

    log_level = logging.INFO
    include_exc_info = False

    def __init__(self, message: str):
        super().__init__(message=message, status_code=status.HTTP_400_BAD_REQUEST)


class VCFDimensionsEntryMissingError(DivBaseAPIException):
    """Raised when there are no entries in the VCF dimensions db table."""

    log_level = logging.INFO
    include_exc_info = False

    def __init__(self, project_name: str | None = None):
        if project_name:
            message = (
                f"The VCF dimensions cache in project '{project_name}' is missing or empty. "
                "Please ensure that there are VCF files in the project and run:\n"
                "'divbase-cli dimensions update --project <project_name>'\n"
            )
        else:
            message = (
                "The VCF dimensions cache is missing or empty. "
                "Please ensure that there are VCF files in the project and run:\n"
                "'divbase-cli dimensions update --project <project_name>'\n"
            )
        super().__init__(message, status_code=status.HTTP_404_NOT_FOUND)


class DimensionsUpdateAlreadyInProcessError(DivBaseAPIException):
    """Raised when a user tries to queue a new VCF dimensions update task for a project that already has a queued or running update task."""

    log_level = logging.INFO
    include_exc_info = False

    def __init__(self, project_name: str, ongoing_task_id: int):
        message = (
            f"A VCF dimensions update task (with id: {ongoing_task_id}) is already in process for the project '{project_name}'. \n"
            f"Only one dimensions update task can be run at a time. \n"
            "The dimensions update task will incorporate files uploaded/modified during runtime, so there is no need to start another update task while this one is in progress."
        )
        super().__init__(message=message, status_code=status.HTTP_409_CONFLICT)


class DownloadedFileChecksumMismatchError(DivBaseAPIException):
    """Raised when a worker downloads a file but the calculated checksum does not match the checksum provided by s3."""

    log_level = logging.ERROR
    include_exc_info = True

    def __init__(self, file_path: Path, calculated_checksum: str, expected_checksum: str):
        message = (
            f"Downloaded file checksum mismatch for file '{file_path.name}'. "
            f"Expected: {expected_checksum}, "
            f"but calculated: {calculated_checksum}."
        )
        super().__init__(message=message, status_code=status.HTTP_409_CONFLICT)


class TooManyObjectsInRequestError(DivBaseAPIException):
    """Raise when e.g. too many files are requested to be downloaded in a single request."""

    log_level = logging.WARNING
    include_exc_info = False

    def __init__(self, message: str = "Too many objects to work on in a single request"):
        super().__init__(message=message, status_code=status.HTTP_400_BAD_REQUEST)


class ObjectDoesNotExistError(DivBaseAPIException):
    """Raised when an S3 object/key does not exist in the project's bucket."""

    log_level = logging.INFO
    include_exc_info = False

    def __init__(self, key: str, bucket_name: str):
        message = f"The file/object '{key}' does not exist in the project's store: '{bucket_name}'. "
        super().__init__(message=message, status_code=status.HTTP_404_NOT_FOUND)


class TSVFileNotFoundInProjectError(DivBaseAPIException):
    """Raised when a TSV doesn't exist in project storage."""

    log_level = logging.INFO
    include_exc_info = False

    def __init__(self, filename: str | None = None, project_name: str | None = None):
        self.filename = filename
        self.project_name = project_name
        if filename and project_name:
            message = (
                f"The sample metadata TSV file '{filename}' "
                f"was not found in your project '{project_name}'.\n"
                "Please check that you have spelled the file name correctly and that the file has been uploaded to the project."
            )
        else:
            message = (
                "A sample metadata TSV file was not found.\n"
                "Please check that you have spelled the file name correctly and that the file has been uploaded to the project."
            )
        super().__init__(message=message, status_code=status.HTTP_404_NOT_FOUND)


class QueueClosedError(DivBaseAPIException):
    """Raised when the queue is closed for new tasks (e.g. before a maintenance window)."""

    log_level = logging.INFO
    include_exc_info = False

    def __init__(self, message: str):
        super().__init__(message=message, status_code=status.HTTP_400_BAD_REQUEST)


class PATLimitExceededError(DivBaseAPIException):
    """Raised when a user has reached the maximum number of active personal access tokens."""

    log_level = logging.INFO
    include_exc_info = False

    def __init__(
        self,
        message: str = "You have reached the maximum of 5 active tokens. Please revoke one before creating a new one.",
    ):
        super().__init__(message=message, status_code=status.HTTP_409_CONFLICT)


class PATDuplicateNameError(DivBaseAPIException):
    """Raised when a user already has an active (not soft-deleted) personal access token with the same name."""

    log_level = logging.INFO
    include_exc_info = False

    def __init__(self, message: str = "You already have an active token with this name."):
        super().__init__(message=message, status_code=status.HTTP_409_CONFLICT)
