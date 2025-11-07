"""
DivBase API custom exceptions.

All exceptions in this module should inherit from DivBaseAPIException, so we can catch for that externally.
"""

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


class BucketVersioningFileAlreadyExistsError(DivBaseAPIException):
    """Raised when trying to create a bucket versioning file that already exists for that project"""

    def __init__(self, bucket_name: str):
        message = (
            f"The bucket versioning file already exists for the bucket: '{bucket_name}'.\n"
            "You can already add a new bucket version to this file using the 'add' command"
        )
        super().__init__(message, status_code=status.HTTP_400_BAD_REQUEST)


class BucketVersionAlreadyExistsError(DivBaseAPIException):
    """Raised when user tries to add new version with version name same as prexisting version."""

    def __init__(self, version_name: str, bucket_name: str):
        message = f"You're trying to add a version: '{version_name}' that already exists in the bucket: '{bucket_name}'"
        super().__init__(message, status_code=status.HTTP_400_BAD_REQUEST)


class BucketVersionNotFoundError(DivBaseAPIException):
    """Raised when the specified bucket version file is not found in the bucket versioning dictionary"""

    def __init__(self, bucket_version: str, bucket_name: str):
        message = f"The version requested: '{bucket_version}' does not exist in the bucket: '{bucket_name}'."
        super().__init__(message, status_code=status.HTTP_404_NOT_FOUND)


class VCFDimensionsEntryMissingError(DivBaseAPIException):
    """Raised when there are no entries in the VCF dimensions db table."""

    def __init__(self, project_name: str):
        message = (
            f"The VCF dimensions index in project '{project_name}' is missing or empty. "
            "Please ensure that there are VCF files in the project and run:\n"
            "'divbase-cli dimensions update --project <project_name>'\n"
        )
        super().__init__(message, status_code=status.HTTP_404_NOT_FOUND)
