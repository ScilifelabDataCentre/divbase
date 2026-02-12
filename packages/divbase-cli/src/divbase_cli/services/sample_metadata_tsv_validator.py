"""
Client-side validator for DivBase sidecar metadata TSV files.

Requires that the dimensions index for the project is up to date and is fetched from the server.
Validates formatting requirements without sending data to the server.
"""

import csv
from pathlib import Path


class MetadataTSVValidator:
    """Validates sidecar metadata TSV files against DivBase requirements."""

    FORBIDDEN_CHARS = [","]

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

        Returns a tuple of (stats, errors, warnings) where stats is a dictionary of collected statistics,
        errors is a list of error messages, and warnings is a list of warning messages.
        """
        validator = cls(file_path, project_samples)

        try:
            with open(validator.file_path, "r", newline="", encoding="utf-8") as f:
                reader = csv.reader(f, delimiter="\t")
                rows = list(reader)
        except Exception as e:
            validator.errors.append(f"Failed to read file: {e}")
            return validator.stats, validator.errors, validator.warnings

        if not rows:
            validator.errors.append("File is empty")
            return validator.stats, validator.errors, validator.warnings

        validator._validate_header(rows[0])

        if len(rows) > 1:
            validator._validate_data_rows(rows)

        return validator.stats, validator.errors, validator.warnings

    def _validate_header(self, header: list[str]) -> None:
        """Validate the header row."""
        if not header:
            self.errors.append("Header row is missing")
            return

        if header[0] != "#Sample_ID":
            self.errors.append(f"First column must be named '#Sample_ID', found: '{header[0]}'")

        if len(header) != len(set(header)):
            duplicates = [col for col in header if header.count(col) > 1]
            self.errors.append(f"Duplicate column names found: {set(duplicates)}")

        for i, col in enumerate(header):
            if not col.strip():
                self.errors.append(f"Empty column name at position {i + 1}")

    def _validate_data_rows(self, rows: list[list[str]]) -> None:
        """Validate all data rows."""
        header = rows[0]
        data_rows = rows[1:]

        num_columns = len(header)
        sample_ids_seen = set()
        tsv_samples = set()

        column_types: dict[int, set[str]] = {i: set() for i in range(1, num_columns)}
        empty_cells_per_column: dict[str, int] = {header[i]: 0 for i in range(1, num_columns)}

        has_multi_values = False

        for row_num, row in enumerate(data_rows, start=2):  # Start at row 2 (after header)
            if len(row) != num_columns:
                sample_hint = f" (Sample_ID: '{row[0]}')" if row else ""
                self.errors.append(
                    f"Row {row_num}: Expected {num_columns} tab-separated columns from reading the header, found {len(row)}{sample_hint}. "
                    "Check that all values are separated by tabs (not spaces)."
                )
                continue

            sample_id = row[0].strip() if row else ""

            if not sample_id:
                self.errors.append(f"Row {row_num}: Sample_ID is empty")
                continue

            if ";" in sample_id:
                self.errors.append(
                    f"Row {row_num}: Sample_ID '{sample_id}' contains semicolon. Sample_ID must contain only one value."
                )

            if sample_id in sample_ids_seen:
                self.errors.append(f"Row {row_num}: Duplicate Sample_ID '{sample_id}'")
            else:
                sample_ids_seen.add(sample_id)
                tsv_samples.add(sample_id)

            for col_idx, cell in enumerate(row):
                self._validate_cell(row_num, col_idx, header[col_idx], cell)

                # Track column types and empty-cells for user-defined columns (skip col 0, i.e. Sample_ID)
                if col_idx > 0:
                    if not cell.strip():
                        empty_cells_per_column[header[col_idx]] += 1

                    if ";" in cell:
                        has_multi_values = True
                    self._infer_column_type(row_num, col_idx, header[col_idx], cell, column_types)

        self._check_mixed_types(header, column_types)

        self._validate_sample_names(tsv_samples)

        self._collect_statistics(header, tsv_samples, column_types, has_multi_values, empty_cells_per_column)

    def _validate_cell(self, row_num: int, col_idx: int, col_name: str, cell: str) -> None:
        """Validate an individual cell."""

        for char in self.FORBIDDEN_CHARS:
            if char in cell:
                self.errors.append(f"Row {row_num}, Column '{col_name}': Cell contains forbidden character '{char}'")

        if cell != cell.strip():
            self.warnings.append(
                f"Row {row_num}, Column '{col_name}': Cell has leading or trailing whitespace "
                "(will be stripped by server)"
            )

        if col_idx > 0 and not cell.strip():
            self.warnings.append(
                f"Row {row_num}, Column '{col_name}': Cell is empty. "
                "Empty values will be treated as missing by the server and will not match any filter conditions in queries."
            )

    def _infer_column_type(
        self, row_num: int, col_idx: int, col_name: str, cell: str, column_types: dict[int, set[str]]
    ) -> None:
        """
        Infer the type of values in a column and track type information.
        Matches server-side logic in queries.py::_is_semicolon_separated_numeric_column

        Columns with a mix of numeric-looking and non-numeric values (e.g., "8", "1a", "5a")
        are treated as string columns. Mixed types are reported as warnings, not errors.
        """
        values = [v.strip() for v in cell.split(";") if v.strip()]

        cell_has_numeric = False
        cell_has_string = False

        for value in values:
            # Try to determine if numeric or string first
            # Note! The queries use Pandas for this which is not used here due to different dependencies in the packages. There could potentially be a discrepancy here.
            try:
                float(value)
                cell_has_numeric = True
                column_types[col_idx].add("numeric")
            except ValueError:
                cell_has_string = True
                column_types[col_idx].add("string")

                # Check for hyphens in non-numeric values that might indicate range notation.
                # Negative numbers should already have been classified as numeric with the float() check.
                # This is a warning to help users who may have used range notation (e.g., "1-2") instead of
                # semicolons (e.g., "1;2") in their data values.
                if (
                    "-" in value
                    and any(c.isdigit() for c in value)
                    and ("numeric" in column_types[col_idx] or all(t == "numeric" for t in column_types[col_idx] if t))
                ):
                    self.warnings.append(
                        f"Row {row_num}, Column '{col_name}': Value '{value}' contains a hyphen. "
                        f"This appears to be range notation (e.g., '1-2'), which is not allowed in data values. "
                        f"If this is meant to be a numeric multi-value column, use semicolons to separate values (e.g., '1;2'). "
                        f"This column will be treated as a string column."
                    )

        # Check for mixed types within the same cell (e.g., "1;abc") and warn the user if applicable
        if cell_has_numeric and cell_has_string:
            self.warnings.append(
                f"Row {row_num}, Column '{col_name}': Cell '{cell}' contains mixed types "
                f"(both numeric and non-numeric values in semicolon-separated cell). "
                f"This column will be treated as a string column."
            )

    def _check_mixed_types(self, header: list[str], column_types: dict[int, set[str]]) -> None:
        """
        Check for mixed types in columns and report as informational warnings.
        Matches server-side logic in queries.py::_is_semicolon_separated_numeric_column

        Columns with mixed types are treated as string columns by the DivBase query engine.
        This happen for values such as e.g., "8", "1a", "5a" that happen to look numeric but
        are semantically a strings (e.g. names, IDs)..
        """
        mixed_columns = []
        for col_idx, types in column_types.items():
            if len(types) > 1:
                col_name = header[col_idx]
                mixed_columns.append(col_name)

        if mixed_columns:
            self.warnings.append(
                f"The following columns contain mixed types (both numeric-looking and string values): {mixed_columns}. "
                "A column is only numeric if all values (including each part in semicolon-separated cells) are valid numbers. "
                "These columns will be treated as string columns by DivBase. Numeric query operations "
                "(ranges, inequalities) will not be applicable to these columns."
            )

    def _validate_sample_names(self, tsv_samples: set[str]) -> None:
        """Validate sample names against project dimensions."""

        missing_from_project = tsv_samples - self.project_samples
        if missing_from_project:
            examples = sorted(list(missing_from_project))
            self.errors.append(
                f"The following samples in the TSV were not found in the DivBase project's dimensions index: {examples}. "
                "DivBase requires that all samples in the TSV file must be present in the project's dimensions index to be used for queries."
            )

        missing_from_tsv = self.project_samples - tsv_samples
        if missing_from_tsv:
            examples = sorted(list(missing_from_tsv))
            self.warnings.append(
                f"The following samples in the DivBase project's dimensions index were not found in the TSV: {examples}. "
                "This is allowed for DivBase metadata TSV files, but please be aware that these samples will not be considered when making queries with this metadata file."
            )

    def _collect_statistics(
        self,
        header: list[str],
        tsv_samples: set[str],
        column_types: dict[int, set[str]],
        has_multi_values: bool,
        empty_cells_per_column: dict[str, int],
    ) -> None:
        """Collect statistics about the TSV file."""

        self.stats["total_columns"] = len(header)
        self.stats["user_defined_columns"] = len(header) - 1  # Exclude Sample_ID

        matching_samples = tsv_samples & self.project_samples
        self.stats["samples_in_tsv"] = len(tsv_samples)
        self.stats["samples_matching_project"] = len(matching_samples)
        self.stats["total_project_samples"] = len(self.project_samples)

        numeric_cols = []
        string_cols = []
        mixed_cols = []

        for col_idx, types in column_types.items():
            col_name = header[col_idx]
            if len(types) > 1:
                mixed_cols.append(col_name)
            elif "numeric" in types:
                numeric_cols.append(col_name)
            elif "string" in types:
                string_cols.append(col_name)

        self.stats["numeric_columns"] = numeric_cols
        self.stats["string_columns"] = string_cols
        self.stats["mixed_type_columns"] = mixed_cols
        self.stats["has_multi_values"] = has_multi_values

        # Only include columns with empty cells in stats
        columns_with_empty_cells = {col: count for col, count in empty_cells_per_column.items() if count > 0}
        if columns_with_empty_cells:
            self.stats["empty_cells_per_column"] = columns_with_empty_cells
