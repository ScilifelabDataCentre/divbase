"""
Custom exceptions for DivBase packages.

These are raised by lover-level functions/methods which understand the context of the error.

Note: By adding the `__str__` method to each exception,
we ensure that when you manually raise a specific exception the error message looks good
"""

from pathlib import Path

from divbase_lib.api_schemas.s3 import ExistingFileResponse


class ObjectDoesNotExistError(FileNotFoundError):
    """Raised when an S3 object/key does not exist in the bucket."""

    def __init__(self, key: str, bucket_name: str):
        error_message = f"The file/object '{key}' does not exist in the bucket '{bucket_name}'. "
        super().__init__(error_message)
        self.key = key
        self.bucket = bucket_name
        self.error_message = error_message

    def __str__(self):
        return self.error_message


class ObjectDoesNotExistInSpecifiedVersionError(KeyError):
    """Raised when an S3 object/key does not exist in the specified bucket versioning yaml file."""

    def __init__(self, project_name: str, bucket_version: str, missing_objects: list[str]):
        missing_objects_str = "\n".join(f"- '{name}'" for name in missing_objects)
        error_message = (
            f"For the project: '{project_name}'\n"
            f"And bucket version you specified: '{bucket_version}':\n"
            "The following objects could not be found in the metadata file:\n"
            f"{missing_objects_str}"
            "\n Maybe they only existed in a later version of the project/bucket?"
        )
        self.bucket_name = project_name
        self.bucket_version = bucket_version
        self.missing_objects = missing_objects
        self.error_message = error_message

    def __str__(self):
        return self.error_message


class BucketVersioningFileDoesNotExist(FileNotFoundError):
    """
    Raised when the bucket versioning file does not exist in the bucket,
    and it needs to be for the given operation (e.g. add/delete version).
    """

    def __init__(self, bucket_name: str):
        error_message = (
            f"The bucket: '{bucket_name}', does not have a bucket versioning file.\n"
            "please create one first using the 'divbase version create' command."
        )
        super().__init__(error_message)
        self.bucket_name = bucket_name
        self.error_message = error_message

    def __str__(self):
        return self.error_message


class FilesAlreadyInBucketError(FileExistsError):
    """
    Raised when trying to upload file(s) that already exists in the bucket
    and the user does not want to accidently create a new version of any file.
    """

    def __init__(self, existing_objects: list[ExistingFileResponse], project_name: str):
        files_list = "\n".join(f"- '{obj.object_name}'" for obj in existing_objects)
        error_message = (
            f"For the project: '{project_name}'\n"
            "The following objects that you're trying to upload already exist in the project's bucket:\n"
            f"{files_list}."
        )
        super().__init__(error_message)
        self.existing_objects = existing_objects
        self.project_name = project_name
        self.error_message = error_message

    def __str__(self):
        return self.error_message


class ProjectNameNotSpecifiedError(Exception):
    """
    Raised when the project name is not specified in the command line arguments, and
    no default project is set in the user config file.
    """

    def __init__(self, config_path: Path):
        error_message = (
            "No project name provided. \n"
            f"Please either set a default project in your user configuration file at '{config_path.resolve()}'.\n"
            f"or pass the flag '--project <project_name>' to this command.\n"
        )
        super().__init__(error_message)
        self.config_path = config_path
        self.error_message = error_message

    def __str__(self):
        return self.error_message


class ProjectNotInConfigError(Exception):
    """
    Raised when the project name was
        1. specified in the command line arguments OR
        2. set as the default project in the user config file.
    But info about the project could not be obtained from the user config file.
    """

    def __init__(self, config_path: Path, project_name: str):
        error_message = (
            f"Couldn't get information about the project named: '{project_name}' \n"
            f"Please check the project is included in '{config_path.resolve()}'.\n"
            f"you can run 'divbase-cli config show' to view the contents of your config file.\n"
        )
        super().__init__(error_message)
        self.config_path = config_path
        self.project_name = project_name
        self.error_message = error_message

    def __str__(self):
        return self.error_message


class BcftoolsEnvironmentError(Exception):
    """Raised when there's an issue with the execution environment (Docker, etc.)."""

    def __init__(self, container_name: str):
        self.container_name = container_name
        error_message = (
            f"No running container found with name {self.container_name}. Ensure the Docker image is available.\n"
        )
        super().__init__(error_message)

        self.error_message = error_message

    def __str__(self):
        return self.error_message


class BcftoolsCommandError(Exception):
    """Raised when a bcftools command fails to execute properly."""

    def __init__(self, command: str, error_details: Exception = None):
        self.command = command
        self.error_details = error_details

        error_message = f"bcftools command failed: '{command}'"
        if error_details:
            error_message += f" with error details: {error_details}"

        super().__init__(error_message)

    def __str__(self):
        if hasattr(self.error_details, "stderr") and self.error_details.stderr:
            return f"bcftools command failed: '{self.command}' with error: {self.error_details.stderr}"
        return super().__str__()


class BcftoolsPipeEmptyCommandError(Exception):
    """Raised when an empty command is provided to the bcftools pipe."""

    def __init__(self):
        error_message = "Empty command provided. Please specify at least one valid bcftools command."
        super().__init__(error_message)
        self.error_message = error_message

    def __str__(self):
        return self.error_message


class BcftoolsPipeUnsupportedCommandError(Exception):
    """Raised when a bcftools command unsupported by the BcftoolsQueryManager class is provided."""

    def __init__(self, command: str, position: int, valid_commands: list[str]):
        self.command = command
        self.position = position
        self.valid_commands = valid_commands

        message = (
            f"Unsupported bcftools command '{command}' at position {position}. "
            f"Only the following commands are supported: {', '.join(valid_commands)}"
        )
        super().__init__(message)


class SidecarNoDataLoadedError(Exception):
    """Raised when no data is loaded in SidecarQueryManager."""

    def __init__(self, file_path: Path, submethod: str, error_details: str | None = None):
        self.file_path = file_path
        self.error_details = error_details

        error_message = f"No data loaded from file '{file_path}', as raised in submethod '{submethod}'."
        if error_details:
            error_message += f"More details about the error: {error_details}"
        super().__init__(error_message)
        self.error_message = error_message

    def __str__(self):
        return self.error_message


class SidecarInvalidFilterError(Exception):
    """Raised when an invalid filter is provided to SidecarQueryManager."""

    pass


class SidecarColumnNotFoundError(Exception):
    """Raised when a requested column is not found in the query result."""

    pass


class NoVCFFilesFoundError(Exception):
    """Raised when no VCF files are found in the project bucket."""

    pass


class ChecksumVerificationError(Exception):
    """Raised when a calculated file's checksum does not match the expected value."""

    def __init__(self, expected_checksum: str, calculated_checksum: str):
        self.expected_checksum = expected_checksum
        self.calculated_checksum = calculated_checksum

        message = f"Checksum verification failed. Expected: {expected_checksum}, Calculated: {calculated_checksum}"
        super().__init__(message)
