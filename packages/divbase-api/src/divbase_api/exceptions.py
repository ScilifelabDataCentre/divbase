"""
DivBase API custom exceptions.

All exceptions in this module should inherit from DivBaseAPIException, so we can catch for that externally.
"""

from pathlib import Path

from fastapi import status


class DivBaseAPIException(Exception):
    """Base exception for all DivBase errors."""

    def __init__(self, message: str, status_code: int, headers: dict[str, str] | None = None):
        self.message = message
        self.status_code = status_code
        self.headers = headers or {}
        super().__init__(self.message)


class AuthenticationError(DivBaseAPIException):
    def __init__(self, message: str = "Authentication required"):
        default_headers = {"WWW-Authenticate": "Bearer"}
        super().__init__(message=message, status_code=status.HTTP_401_UNAUTHORIZED, headers=default_headers)


class AuthorizationError(DivBaseAPIException):
    def __init__(self, message: str = "Authorization required"):
        default_headers = {"WWW-Authenticate": "Bearer"}
        super().__init__(message=message, status_code=status.HTTP_403_FORBIDDEN, headers=default_headers)


class UserRegistrationError(DivBaseAPIException):
    """Exception for failed registration of new user or update of existing user."""

    def __init__(
        self,
        internal_logging_message: str,  # for logs, can be more detailed
        user_message: str = "Registration failed. Please try again.",  # given to end user, generic
    ):
        self.user_message = user_message
        super().__init__(message=internal_logging_message, status_code=status.HTTP_400_BAD_REQUEST)


class ProjectNotFoundError(DivBaseAPIException):
    def __init__(self, message: str = "Project not found or you don't have access"):
        super().__init__(message=message, status_code=status.HTTP_404_NOT_FOUND)


class ProjectMemberNotFoundError(DivBaseAPIException):
    def __init__(self, message: str = "Project member not found"):
        super().__init__(message=message, status_code=status.HTTP_404_NOT_FOUND)


class ProjectCreationError(DivBaseAPIException):
    def __init__(self, message: str = "Project creation failed"):
        super().__init__(message=message, status_code=status.HTTP_400_BAD_REQUEST)


class TooManyObjectsInRequestError(DivBaseAPIException):
    """Raise when e.g. too many files are requested to be downloaded in a single request."""

    def __init__(self, message: str = "Too many objects to work on in a single request"):
        super().__init__(message=message, status_code=status.HTTP_400_BAD_REQUEST)


class ProjectVersionCreationError(DivBaseAPIException):
    """Raised when there is an error creating a new project version."""

    def __init__(self, message: str = "Failed to create a new project version"):
        super().__init__(message=message, status_code=status.HTTP_400_BAD_REQUEST)


class ProjectVersionAlreadyExistsError(DivBaseAPIException):
    """Raised when attempting to create a project version that already exists."""

    def __init__(
        self, message: str = "A project version with the specified name already exists, please choose a different name."
    ):
        super().__init__(message=message, status_code=status.HTTP_400_BAD_REQUEST)


class ProjectVersionNotFoundError(DivBaseAPIException):
    """Raised when a project version is not found."""

    def __init__(self, message: str = "Could not find the specified project version"):
        super().__init__(message=message, status_code=status.HTTP_404_NOT_FOUND)


class VCFDimensionsEntryMissingError(DivBaseAPIException):
    """Raised when there are no entries in the VCF dimensions db table."""

    def __init__(self, project_name: str):
        message = (
            f"The VCF dimensions index in project '{project_name}' is missing or empty. "
            "Please ensure that there are VCF files in the project and run:\n"
            "'divbase-cli dimensions update --project <project_name>'\n"
        )
        super().__init__(message, status_code=status.HTTP_404_NOT_FOUND)


class TaskNotFoundInBackendError(DivBaseAPIException):
    """Raised when a task ID exists in the task_history table in the database but not in the results backend."""

    def __init__(
        self,
        message: str = "Task ID not found in results backend. It may have been purged during cleanup of old task records.",
    ):
        super().__init__(message=message, status_code=status.HTTP_410_GONE)


class DownloadedFileChecksumMismatchError(DivBaseAPIException):
    """Raised when a worker downloads a file but the calculated checksum does not match the checksum provided by s3."""

    def __init__(self, file_path: Path, calculated_checksum: str, expected_checksum: str):
        message = (
            f"Downloaded file checksum mismatch for file '{file_path.name}'. "
            f"Expected: {expected_checksum}, "
            f"but calculated: {calculated_checksum}."
        )
        super().__init__(message=message, status_code=status.HTTP_409_CONFLICT)
