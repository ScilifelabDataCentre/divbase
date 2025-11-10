"""
Custom exceptions for the divbase CLI.
"""


class DivBaseCLIError(Exception):
    """Base exception for all divbase CLI errors."""

    pass


class AuthenticationError(DivBaseCLIError):
    """Raised for user authentication errors when using CLI tool."""

    def __init__(self, error_message: str = "Authentication required, make sure you're logged in."):
        super().__init__(error_message)


class DivBaseAPIConnectionError(DivBaseCLIError):
    """Raised when CLI tool can't connect to the provided DivBase API URL."""

    def __init__(
        self,
        error_message: str = "Unable to connect to the DivBase API. Check the API URL and your network connection. Perhaps the server is down?",
    ):
        super().__init__(error_message)


class DivBaseAPIError(DivBaseCLIError):
    """
    Used by CLI tool when making requests to DivBase API.
    Raised when the DivBase API/server responds with an error status code.
    Provides a helpful and easy-to-read error message for the user.
    """

    def __init__(
        self,
        error_details: str = "Not Provided",
        error_type: str = "unknown",
        status_code: int = None,
        http_method: str = "unknown",
        url: str = "unknown",
    ):
        self.status_code = status_code
        self.error_type = error_type
        self.error_details = error_details
        self.http_method = http_method
        self.url = url

        self.error_message = (
            f"DivBase Server returned an error response:\n"
            f"HTTP Status code: {status_code}\n"
            f"HTTP method: {http_method}\n"
            f"URL: {url}\n"
            f"Error type: {error_type}\n"
            f"Details: {error_details}\n"
        )
        super().__init__(self.error_message)

    def __str__(self):
        return self.error_message
