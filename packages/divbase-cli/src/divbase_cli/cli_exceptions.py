"""
Custom exceptions for the divbase CLI.
"""

from pathlib import Path

from divbase_lib.api_schemas.s3 import ExistingFileResponse


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
        error_message = (
            f"DivBase Server returned an error response:\n"
            f"HTTP Status code: {status_code}\n"
            f"HTTP method: {http_method}\n"
            f"URL: {url}\n"
            f"Error type: {error_type}\n"
            f"Details: {error_details}\n"
        )
        self.error_message = error_message
        super().__init__(error_message)


class FileDoesNotExistInSpecifiedVersionError(DivBaseCLIError):
    """Raised when a file does not exist in the project at the specified project version"""

    def __init__(self, project_name: str, project_version: str, missing_files: list[str]):
        missing_files_str = "\n".join(f"- '{name}'" for name in missing_files)
        self.project_name = project_name
        self.project_version = project_version
        self.missing_files = missing_files

        error_message = (
            f"For the project: '{project_name}'\n"
            f"And project version you specified: '{project_version}':\n"
            "The following file(s) could not be found:\n"
            f"{missing_files_str}"
            "\n Maybe they only existed in a later version of the project?"
        )
        super().__init__(error_message)


class FilesAlreadyInProjectError(DivBaseCLIError):
    """
    Raised when trying to upload file(s) that already exists in the project
    and the user does not want to accidently create a new version of any file.
    """

    def __init__(self, existing_files: list[ExistingFileResponse], project_name: str):
        files_list = "\n".join(f"- '{obj.object_name}'" for obj in existing_files)
        self.existing_files = existing_files
        self.project_name = project_name

        error_message = (
            f"For the project: '{project_name}'\n"
            "The exact version of the following file(s) that you're trying to upload already exist inside the project:\n"
            f"{files_list}."
        )
        super().__init__(error_message)


class ProjectNameNotSpecifiedError(DivBaseCLIError):
    """
    Raised when the project name is not specified in the command line arguments, and
    no default project is set in the user config file.
    """

    def __init__(self, config_path: Path):
        self.config_path = config_path
        error_message = (
            "No project name provided. \n"
            f"Please either set a default project in your user configuration file at '{config_path.resolve()}'.\n"
            f"or pass the flag '--project <project_name>' to this command.\n"
        )
        super().__init__(error_message)


class ProjectNotInConfigError(DivBaseCLIError):
    """
    Raised when the project name was
        1. specified in the command line arguments OR
        2. set as the default project in the user config file.
    But info about the project could not be obtained from the user config file.
    """

    def __init__(self, config_path: Path, project_name: str):
        self.config_path = config_path
        self.project_name = project_name
        error_message = (
            f"Couldn't get information about the project named: '{project_name}' \n"
            f"Please check the project is included in '{config_path.resolve()}'.\n"
            f"you can run 'divbase-cli config show' to view the contents of your config file.\n"
        )
        super().__init__(error_message)


class ConfigFileNotFoundError(DivBaseCLIError):
    """Raised when the user's config file cannot be found."""

    def __init__(
        self,
        error_message: str = (
            "You're DivBase configuration file was not found or does not exist.\n"
            "To create a user configuration file, run 'divbase-cli config create'.\n"
            "If you already have a user configuration file that but it is not stored in the default location, "
            "you can pass the '--config <path>' flag to specify the location. \n"
            "You very probably want to just run 'divbase-cli config create' though."
        ),
    ):
        super().__init__(error_message)
