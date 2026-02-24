"""
Client-side validator for DivBase sidecar metadata TSV files.

Requires that the dimensions index for the project is up to date and is fetched from the server.
Validates formatting requirements without sending data to the server.

Calls the shared sidecar metadata TSV validation logic from divbase_lib.metadata_validator.
"""

from pathlib import Path

from divbase_lib.metadata_validator import SharedMetadataValidator


class MetadataTSVValidator:
    """Validates sidecar metadata TSV files against DivBase requirements."""

    def __init__(self, file_path: Path, project_samples: list[str] | set[str]):
        """
        Initialize the validator. File path is the path to the TSV file to validate,
        and project_samples is a list or set of unique sample names from the project's dimensions index.
        """
        self.file_path = file_path
        self.project_samples = set(project_samples) if isinstance(project_samples, list) else project_samples
        self.errors: list[str] = []
        self.warnings: list[str] = []
        self.stats: dict = {}

    @classmethod
    def validate(cls, file_path: Path, project_samples: list[str] | set[str]) -> tuple[dict, list[str], list[str]]:
        """
        Validate a TSV file and return results.

        Since this is run on the client-side, it should try to detect all formatting errors and warnings.
        Therefore, assuming that it is able read the local file, it does not raise exceptions like the SidecarQueryManager does for server-side validation,
        but instead continues working through the file to collects all encountered errors and warnings.

        Returns a tuple of (stats, errors, warnings) where stats is a dictionary of collected statistics,
        errors is a list of error messages, and warnings is a list of warning messages.
        """
        validator = cls(file_path, project_samples)

        shared_validator = SharedMetadataValidator(
            file_path=file_path,
            project_samples=validator.project_samples,
            skip_dimensions_check=False,
        )
        result = shared_validator.load_and_validate()

        validator.errors = result.errors
        validator.warnings = result.warnings

        if result.df is not None and "Sample_ID" in result.df.columns:
            try:
                tsv_samples = set(result.df["Sample_ID"].tolist())
                validator._collect_statistics(
                    df=result.df,
                    tsv_samples=tsv_samples,
                    numeric_cols=result.numeric_columns,
                    string_cols=result.string_columns,
                    mixed_cols=result.mixed_type_columns,
                )
            except (AttributeError, TypeError):
                # If Sample_ID access fails (e.g., in the very rare case that duplicate Sample_ID column make it a DataFrame due to dataframe nesting), skip statistics as the validation errors already captured the issue
                pass

        return validator.stats, validator.errors, validator.warnings

    def _collect_statistics(
        self,
        df,
        tsv_samples: set[str],
        numeric_cols: list[str],
        string_cols: list[str],
        mixed_cols: list[str],
    ) -> None:
        """Collect statistics about the TSV file."""

        self.stats["total_columns"] = len(df.columns)
        self.stats["user_defined_columns"] = len(df.columns) - 1  # Exclude Sample_ID

        matching_samples = tsv_samples & self.project_samples
        self.stats["samples_in_tsv"] = len(tsv_samples)
        self.stats["samples_matching_project"] = len(matching_samples)
        self.stats["total_project_samples"] = len(self.project_samples)

        self.stats["numeric_columns"] = numeric_cols
        self.stats["string_columns"] = string_cols
        self.stats["mixed_type_columns"] = mixed_cols

        has_multi_values = False
        for col in df.columns:
            for val in df[col].dropna():
                if isinstance(val, list):
                    has_multi_values = True
                    break
            if has_multi_values:
                break
        # If has_multi_values is True: at least one cell in the DataFrame contains a Python list (multi-value cell).
        self.stats["has_multi_values"] = has_multi_values

        empty_cells_per_column = {}
        for col in df.columns:
            if col != "Sample_ID":
                empty_count = df[col].isna().sum() + (df[col] == "").sum()
                if empty_count > 0:
                    empty_cells_per_column[col] = int(empty_count)

        if empty_cells_per_column:
            self.stats["empty_cells_per_column"] = empty_cells_per_column
