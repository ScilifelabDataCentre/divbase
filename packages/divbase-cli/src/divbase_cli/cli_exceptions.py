"""
Custom exceptions for the divbase CLI.
"""


class DivBaseCLIError(Exception):
    """Base exception for divbase CLI errors."""

    pass


class ChecksumVerificationError(DivBaseCLIError):
    """Raised when a calculated file's checksum does not match the expected value."""

    def __init__(self, expected_checksum: str, calculated_checksum: str):
        self.expected_checksum = expected_checksum
        self.calculated_checksum = calculated_checksum

        message = f"Checksum verification failed. Expected: {expected_checksum}, Calculated: {calculated_checksum}"
        super().__init__(message)
