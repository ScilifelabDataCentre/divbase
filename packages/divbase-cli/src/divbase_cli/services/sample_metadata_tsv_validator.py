"""
Client-side validator for DivBase sidecar metadata TSV files.

Requires that the dimensions index for the project is up to date and is fetched from the server.
Validates formatting requirements without sending data to the server.

Calls the shared sidecar metadata TSV validation logic from divbase_lib.metadata_validator.
"""

from pathlib import Path

from divbase_lib.metadata_validator import SharedMetadataValidator


class ClientSideMetadataTSVValidator:
    """
    Client-side wrapper for validating sidecar metadata TSV files against DivBase requirements.
    Calls SharedMetadataValidator to perform the actual validation.
    """

    def __init__(
        self,
        file_path: Path,
        project_samples: list[str] | set[str],
        dimensions_sample_preview_limit: int | None = 20,
    ):
        """
        Initialize the validator. File path is the path to the TSV file to validate,
        and project_samples is a list or set of unique sample names from the project's dimensions index.
        """
        self.file_path = file_path
        self.project_samples = set(project_samples) if isinstance(project_samples, list) else project_samples
        self.dimensions_sample_preview_limit = dimensions_sample_preview_limit
        self.errors: list[str] = []
        self.warnings: list[str] = []
        self.stats: dict = {}

    def validate(self) -> tuple[dict, list[str], list[str]]:
        """
        Validate a TSV file and return results.

        Since this is run on the client-side, it should try to detect all formatting errors and warnings.
        Therefore, assuming that it is able read the local file, it does not raise exceptions like the SidecarQueryManager does for server-side validation,
        but instead continues working through the file to collects all encountered errors and warnings.

        Returns a tuple of (stats, errors, warnings) where stats is a dictionary of collected statistics,
        errors is a list of error messages, and warnings is a list of warning messages.
        """
        shared_validator = SharedMetadataValidator(
            file_path=self.file_path,
            project_samples=self.project_samples,
            skip_dimensions_check=False,
            dimensions_sample_preview_limit=self.dimensions_sample_preview_limit,
        )
        result = shared_validator.load_and_validate()

        self.errors = [error_entry.message for error_entry in result.errors]
        self.warnings = [warning_entry.message for warning_entry in result.warnings]
        self.stats = result.stats

        return self.stats, self.errors, self.warnings
