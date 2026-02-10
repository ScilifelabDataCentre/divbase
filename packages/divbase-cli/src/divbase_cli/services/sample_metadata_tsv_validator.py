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

    def __init__(self, file_path: Path, project_samples: set[str]):
        """
        Initialize the validator. File path is the path to the TSV file to validate,
        and project_samples is a set of unique sample names from the project's dimensions index.
        """
        self.file_path = file_path
        self.project_samples = project_samples
        self.errors: list[str] = []
        self.warnings: list[str] = []
        self.stats: dict = {}

    def validate(self) -> tuple[dict, list[str], list[str]]:
        """
        Run all validation checks on the TSV file.
        Returns a tuple of (stats, errors, warnings) where stats is a dictionary of collected statistics about the TSV file,
        errors is a list of error messages, and warnings is a list of warning messages.
        """
        self.errors = []
        self.warnings = []
        self.stats = {}

        try:
            with open(self.file_path, "r", newline="", encoding="utf-8") as f:
                reader = csv.reader(f, delimiter="\t")
                rows = list(reader)
        except Exception as e:
            self.errors.append(f"Failed to read file: {e}")
            return self.stats, self.errors, self.warnings

        if not rows:
            self.errors.append("File is empty")
            return self.stats, self.errors, self.warnings

        self._validate_header(rows[0])

        if len(rows) > 1:
            self._validate_data_rows(rows)

        return self.stats, self.errors, self.warnings

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

                # Track column types for user-defined columns (skip col 0, i.e. Sample_ID)
                if col_idx > 0:
                    if ";" in cell:
                        has_multi_values = True
                    self._infer_column_type(row_num, col_idx, header[col_idx], cell, column_types)

        self._check_mixed_types(header, column_types)

        self._validate_sample_names(tsv_samples)

        self._collect_statistics(header, tsv_samples, column_types, has_multi_values)

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

    def _infer_column_type(
        self, row_num: int, col_idx: int, col_name: str, cell: str, column_types: dict[int, set[str]]
    ) -> None:
        """
        Infer the type of values in a column and validate type consistency.
        Matches server-side logic in queries.py::_is_semicolon_separated_numeric_column
        """
        values = [v.strip() for v in cell.split(";") if v.strip()]

        cell_has_numeric = False
        cell_has_string = False

        for value in values:
            # Check for hyphens in values that might be numeric
            # The server-side logic rejects hyphens in numeric columns (e.g., "1-2" could be confused with range syntax)
            if (
                "-" in value
                and any(c.isdigit() for c in value)
                and ("numeric" in column_types[col_idx] or all(t == "numeric" for t in column_types[col_idx] if t))
            ):
                self.errors.append(
                    f"Row {row_num}, Column '{col_name}': Value '{value}' contains a hyphen. "
                    f"Hyphens are not allowed in numeric column values (only in string columns). "
                    f"If this is meant to be a string column, all values should be non-numeric strings."
                )

            # Try to determine if numeric or string. Note! The queries used Pandas for this, so there could potentially be a discrepency here.
            try:
                float(value)
                cell_has_numeric = True
                column_types[col_idx].add("numeric")
            except ValueError:
                cell_has_string = True
                column_types[col_idx].add("string")

        # Check for mixed types within the same cell (e.g., "1;abc")
        if cell_has_numeric and cell_has_string:
            self.errors.append(
                f"Row {row_num}, Column '{col_name}': Cell '{cell}' contains mixed types. "
                f"All cell values in the same column must be consistently numeric or string."
            )

    def _check_mixed_types(self, header: list[str], column_types: dict[int, set[str]]) -> None:
        """
        Check for mixed types in columns and raise errors.
        Matches server-side logic in queries.py::_is_semicolon_separated_numeric_column
        """
        for col_idx, types in column_types.items():
            if len(types) > 1:
                col_name = header[col_idx]
                self.errors.append(
                    f"Column '{col_name}': Contains mixed types (both numeric and string values). "
                    f"All values in a column must be consistently numeric or string for DivBase sidecar metadata queries to work correctly."
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
        self, header: list[str], tsv_samples: set[str], column_types: dict[int, set[str]], has_multi_values: bool
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
