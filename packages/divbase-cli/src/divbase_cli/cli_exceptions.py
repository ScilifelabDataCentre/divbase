"""
Custom exceptions for the divbase CLI.
"""

from pathlib import Path

from divbase_lib.divbase_constants import SUPPORTED_DIVBASE_FILE_TYPES, UNSUPPORTED_CHARACTERS_IN_FILENAMES


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
        status_code: int = 500,
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

    def __init__(self, existing_files: dict[Path, str], project_name: str):
        files_list = "\n".join(
            f"'{file_path}' (Checksum: {checksum})" for file_path, checksum in existing_files.items()
        )
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

    def __init__(self):
        error_message = (
            "No project name provided.\n"
            "Please either set a default project in your user configuration file.\n"
            "or pass the flag '--project <project_name>' to this command.\n"
            "To set a default project, you can run 'divbase-cli config set-default <project_name>'.\n"
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


class UnsupportedFileTypeError(DivBaseCLIError):
    """Raised when one or more files to be uploaded are not supported by DivBase (based on file extension)."""

    def __init__(self, unsupported_files: list[Path], supported_types: tuple[str, ...] = SUPPORTED_DIVBASE_FILE_TYPES):
        self.unsupported_files = unsupported_files
        self.supported_types = supported_types
        message = (
            f"The following file(s) have types that are not supported by DivBase and therefore cannot be uploaded: \n"
            f"{'\n'.join(str(file) for file in unsupported_files)}\n"
            f"DivBase currently supports the following file types: {', '.join(SUPPORTED_DIVBASE_FILE_TYPES)}\n"
            "If you want us to support another file type, please let us know."
        )
        super().__init__(message)


class UnsupportedFileNameError(DivBaseCLIError):
    """Raised when one or more files to be uploaded have unsupported characters in their filenames."""

    def __init__(self, unsupported_files: list[Path]):
        self.unsupported_files = unsupported_files
        message = (
            f"The following file(s) have unsupported characters in their filenames and therefore cannot be uploaded: \n"
            f"{'\n'.join(str(file) for file in unsupported_files)}\n"
            f"Filenames cannot contain any of the following characters: {', '.join(UNSUPPORTED_CHARACTERS_IN_FILENAMES)}\n"
            "Please rename the files and try again."
        )
        super().__init__(message)
