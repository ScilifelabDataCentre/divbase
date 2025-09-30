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
    def __init__(
        self,
        message: str,
        user_message: str = "Registration failed. Please try again.",  # given to end user, generic
    ):
        self.user_message = user_message
        super().__init__(message=message, status_code=status.HTTP_400_BAD_REQUEST)


class ProjectNotFoundError(DivBaseAPIException):
    def __init__(self, message: str = "Project not found or you don't have access"):
        super().__init__(message=message, status_code=status.HTTP_404_NOT_FOUND)


class ProjectMemberNotFoundError(DivBaseAPIException):
    def __init__(self, message: str = "Project member not found"):
        super().__init__(message=message, status_code=status.HTTP_404_NOT_FOUND)


class ProjectCreationError(DivBaseAPIException):
    def __init__(self, message: str = "Project creation failed"):
        super().__init__(message=message, status_code=status.HTTP_400_BAD_REQUEST)
