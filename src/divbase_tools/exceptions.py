"""
Custom exceptions for the divbase_tools package.
"""


class ObjectDoesNotExistError(FileNotFoundError):
    """Raised when an S3 object/key does not exist in the bucket."""

    def __init__(self, key: str, bucket_name: str):
        super().__init__(f"The file '{key}' does not exist in the bucket '{bucket_name}'.")
        self.key = key
        self.bucket = bucket_name
