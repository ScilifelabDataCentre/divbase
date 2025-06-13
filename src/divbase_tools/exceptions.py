"""
Custom exceptions for the divbase_tools package.

These are raised by lover-level functions/methods which understand the context of the error.

Note: By adding the `__str__` method to each exception,
we ensure that when you manually raise a specific exception the error message looks good
"""

from pathlib import Path


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

    def __init__(self, bucket_name: str, bucket_version: str, missing_objects: list[str]):
        missing_objects_str = "\n".join(f"- '{name}'" for name in missing_objects)
        error_message = (
            f"In the bucket: '{bucket_name}'\n"
            f"For the bucket version you specified: '{bucket_version}':\n"
            "The following objects could not be found in the metadata file:\n"
            f"{missing_objects_str}"
            "\n Maybe they only existed in a later version of the bucket?"
        )
        self.bucket_name = bucket_name
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


class BucketVersionNotFoundError(KeyError):
    """Raised when the specified bucket version file is not found in the bucket versioning dictionary"""

    def __init__(self, bucket_version: str, bucket_name: str):
        error_message = (
            f"The version of the bucket specified: '{bucket_version}' does not exist in the bucket '{bucket_name}'."
        )
        super().__init__(error_message)
        self.bucket_version = bucket_version
        self.bucket_name = bucket_name
        self.error_message = error_message

    def __str__(self):
        return self.error_message


class FilesAlreadyInBucketError(FileExistsError):
    """
    Raised when trying to upload file(s) that already exists in the bucket
    and the user does not want to accidently create a new version of any file.

    TODO - This needs some thought in the future, as currenly the path of the file is not used in setting the name of the object,
    only the file name is.
    This decision was taken as the s3 bucket does not have a directory structure.

    But error will be raised if previously uploaded file looks like this: dir1/file1.txt
    and to be uploaded file looks like this: dir2/file1.txt
    """

    def __init__(self, existing_objects: list[str], bucket_name: str):
        files_list = "\n".join(f"- '{name}'" for name in existing_objects)
        error_message = (
            f"In the bucket: '{bucket_name}'\n"
            "The following objects that you're trying to upload already exist in the bucket:\n"
            f"'{files_list}'."
        )
        super().__init__(error_message)
        self.existing_objects = existing_objects
        self.bucket = bucket_name
        self.error_message = error_message

    def __str__(self):
        return self.error_message


class DivBaseCredentialsNotFoundError(Exception):
    """
    Raised when the credentials for accessing DivBase API, which are stored as enviroment variables are not found.
    """

    def __init__(self, access_key_name: str, secret_key_name: str, config_path: Path):
        error_message = (
            "\n"
            "Either your DivBase Access and/or secret key was not found in your enviroment variables'.\n"
            f"Please ensure that the environment variables '{access_key_name}' and '{secret_key_name}' are set.\n"
            f"The simplest way to do this is to create a .env file in the directory you run your commands from.\n"
            f"If you want to use a different name for the environment variables, you can set them in your user config file: {config_path}\n"
        )
        super().__init__(error_message)
        self.config_path = config_path
        self.access_key_name = access_key_name
        self.secret_key_name = secret_key_name
        self.error_message = error_message

    def __str__(self):
        return self.error_message


class BucketNameNotSpecifiedError(Exception):
    """
    Raised when the bucket name is not specified in the command line arguments, and
    no default bucket is set in the user config file.
    """

    def __init__(self, config_path: Path):
        error_message = (
            "No bucket name provided. \n"
            f"Please either set a default bucket in your user configuration file at '{config_path.resolve()}'.\n"
            f"or pass the flag '--bucket-name <bucket_name>' to this command.\n"
        )
        super().__init__(error_message)
        self.config_path = config_path
        self.error_message = error_message

    def __str__(self):
        return self.error_message
