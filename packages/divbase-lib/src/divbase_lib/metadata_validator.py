"""
Shared validation logic for DivBase sidecar metadata TSV files.

This file contains the single source of truth for the TSV content validation logic used by both
the CLI validator (MetadataTSVValidator) on the client-side, and the SidecarQueryManager on the server-side.

Note! Logic for the queries themselves (e.g. how filtering is handled) is not shared between the two.
This file is only for validation of the contents of the TSV file, not for query processing.
"""

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd


@dataclass
class MetadataValidationResult:
    """
    Dataclass to hold the results of the TSV file validation.
    """

    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    stats: dict[str, Any] = field(default_factory=dict)
    df: pd.DataFrame | None = None
    mixed_type_columns: list[str] = field(default_factory=list)
    numeric_columns: list[str] = field(default_factory=list)
    string_columns: list[str] = field(default_factory=list)


class SharedMetadataValidator:
    """
    Core validation logic for DivBase sidecar metadata TSV files. Shared between client-side MetadataTSVValidator
    and server-side SidecarQueryManager to ensure consistent validation behavior.

    It does not validate metadata query filters. That is handled in the SidecarQueryManager.

    This class handles the following validation of the TSV content:
    - Header (duplicates, empty columns, first column name)
    - Sample_ID (empty, duplicates, semicolons)
    - Column type (numeric vs string, semicolon-separated values)
    - Data format (commas, whitespace, column count)
    - Mixed type detection and warnings
    """

    def __init__(self, file_path: Path, project_samples: set[str] | None = None, skip_dimensions_check: bool = False):
        self.file_path = file_path
        self.project_samples = project_samples
        self.skip_dimensions_check = skip_dimensions_check
        self.result = MetadataValidationResult()
        self.df: pd.DataFrame | None = None

    def load_and_validate(self) -> MetadataValidationResult:
        """
        Main entry point to the class. Load a TSV file and call helper methods to validate it.
        """
        try:
            with open(self.file_path, "r", newline="", encoding="utf-8") as f:
                reader = csv.reader(f, delimiter="\t")
                rows = list(reader)

            if not rows:
                self.result.errors.append("File is empty")
                return self.result

            # Pre-pandas checks
            first_line = "\t".join(rows[0])
            header_errors, header_warnings = self._validate_raw_header(first_line)
            self.result.errors.extend(header_errors)
            self.result.warnings.extend(header_warnings)

            if len(rows) > 1:
                row_errors, row_warnings = self._check_row_formatting(rows)
                self.result.errors.extend(row_errors)
                self.result.warnings.extend(row_warnings)

            # Initiate Pandas dataframe from TSV
            df = pd.read_csv(self.file_path, sep="\t", skipinitialspace=True, on_bad_lines="skip")
            df.columns = df.columns.str.lstrip("#")

            for col in df.columns:
                df[col] = df[col].apply(lambda x: x.strip() if isinstance(x, str) else x)

            self.result.df = df
            self.df = df

            # Dataframe checks
            sample_id_errors, sample_id_warnings = self._validate_sample_ids(df)
            self.result.errors.extend(sample_id_errors)
            self.result.warnings.extend(sample_id_warnings)

            mixed_type_columns, cell_warnings = self._detect_mixed_type_columns(df)
            self.result.warnings.extend(cell_warnings)
            self.result.mixed_type_columns = mixed_type_columns

            mixed_type_warning = self._generate_mixed_type_warning(mixed_type_columns)
            if mixed_type_warning:
                self.result.warnings.append(mixed_type_warning)

            comma_warnings = self._check_for_commas(df)
            self.result.warnings.extend(comma_warnings)

            array_notation_warnings = self._check_for_array_notation(df)
            self.result.warnings.extend(array_notation_warnings)

            numeric_cols, string_cols, mixed_cols = self._classify_columns(df, mixed_type_columns)
            self.result.numeric_columns = numeric_cols
            self.result.string_columns = string_cols

            if not self.skip_dimensions_check and self.project_samples is not None:
                tsv_samples = set(df["Sample_ID"].tolist()) if "Sample_ID" in df.columns else set()
                dim_errors, dim_warnings = self._validate_dimensions_match(tsv_samples, self.project_samples)
                self.result.errors.extend(dim_errors)
                self.result.warnings.extend(dim_warnings)

        except Exception as e:
            self.result.errors.append(f"Failed to read file: {e}")

        return self.result

    def _validate_raw_header(self, header_line: str) -> tuple[list[str], list[str]]:
        """
        Validate the raw header line before pandas processing.

        Attempts to catch issues that pandas would silently fix (like duplicate columns).

        """
        errors = []

        raw_columns = header_line.split("\t")
        cleaned_columns = [col.lstrip("#") for col in raw_columns]

        empty_columns = [i + 1 for i, col in enumerate(cleaned_columns) if not col.strip()]
        if empty_columns:
            errors.append(
                f"Empty column name(s) found at position(s): {empty_columns}. All columns must have a non-empty name."
            )

        # Check for duplicate columns after stripping '#'
        seen = {}
        duplicate_columns = []
        for col in cleaned_columns:
            col_stripped = col.strip()
            if col_stripped in seen:
                if col_stripped not in duplicate_columns:
                    duplicate_columns.append(col_stripped)
            else:
                seen[col_stripped] = True

        if duplicate_columns:
            errors.append(
                f"Duplicate column names found: {duplicate_columns}. "
                "Each column name must be unique in the metadata file."
            )

        if raw_columns and raw_columns[0] != "#Sample_ID":
            errors.append(f"First column must be named '#Sample_ID', found: '{raw_columns[0]}'")

        return errors, []

    def _check_row_formatting(self, rows: list[list[str]]) -> tuple[list[str], list[str]]:
        """
        Check for row-level formatting issues that pandas might handle silently.
        """
        errors = []
        warnings = []

        header = rows[0]
        num_columns = len(header)

        for row_num, row in enumerate(rows[1:], start=2):  # Start at row 2 (i.e. skip header)
            if len(row) != num_columns:
                sample_hint = f" (Sample_ID: '{row[0]}')" if row else ""
                errors.append(
                    f"Row {row_num}: Expected {num_columns} tab-separated columns from reading the header, "
                    f"found {len(row)}{sample_hint}. "
                    "Check that all cells in the TSV are separated by tabs (not spaces)."
                )
                continue

            for col_idx, cell in enumerate(row):
                if cell != cell.strip():
                    col_name = header[col_idx]
                    warnings.append(
                        f"Row {row_num}, Column '{col_name}': Cell has leading or trailing whitespace "
                        "(this is allowed, but note that they will be stripped by DivBase server when the TSV is used for queries)"
                    )

                if col_idx == 0 and not cell.strip():
                    errors.append(f"Row {row_num}: Sample_ID is empty")

        return errors, warnings

    def _validate_sample_ids(self, df: pd.DataFrame) -> tuple[list[str], list[str]]:
        """
        Validate Sample_ID column in the DataFrame.

        """
        errors = []

        if "Sample_ID" not in df.columns:
            errors.append("The 'Sample_ID' column is required in the metadata file.")
            return errors, []

        if df["Sample_ID"].isna().any() or (df["Sample_ID"] == "").any():
            errors.append("Sample_ID column contains empty or missing values. All rows must have a valid Sample_ID.")

        if df["Sample_ID"].duplicated().any():
            duplicates = df[df["Sample_ID"].duplicated()]["Sample_ID"].tolist()
            errors.append(f"Duplicate Sample_IDs found: {duplicates}. Each Sample_ID must be unique.")

        semicolon_samples = df[df["Sample_ID"].str.contains(";", na=False)]["Sample_ID"].tolist()
        if semicolon_samples:
            errors.append(
                f"Sample_ID column contains semicolons in values: {semicolon_samples}. "
                "Sample_ID must contain only one value per row (semicolons are not allowed)."
            )

        return errors, []

    def is_semicolon_separated_numeric_column(self, series: pd.Series) -> bool:
        """
        Determine if a column contains semicolon-separated numeric values.

        Pandas infers "1;2;3" as string object dtype. This method checks if all
        non-null values in the column can be parsed as numeric after splitting by semicolon.
        """
        non_null_values = series.dropna()
        if len(non_null_values) == 0:
            return False

        for cell_value in non_null_values:
            if not isinstance(cell_value, str):
                try:
                    float(cell_value)
                    continue
                except (ValueError, TypeError):
                    return False

            # Also check each semicolon-separated part of the cell value
            parts = [p.strip() for p in str(cell_value).split(";") if p.strip()]
            for part in parts:
                try:
                    float(part)
                except ValueError:
                    return False

        return True

    def _detect_mixed_type_columns(self, df: pd.DataFrame) -> tuple[list[str], list[str]]:
        """
        Detect columns with mixed types (numeric and non-numeric values).
        """
        mixed_type_columns = []
        cell_warnings = []

        for col in df.columns:
            if col == "Sample_ID":
                continue

            series = df[col]
            non_null_values = series.dropna()

            if len(non_null_values) == 0:
                continue

            has_numeric = False
            has_string = False

            for idx, cell_value in non_null_values.items():
                if isinstance(cell_value, str) and ";" in cell_value:
                    # Mutli-value cell
                    parts = [p.strip() for p in cell_value.split(";") if p.strip()]
                    cell_has_numeric = False
                    cell_has_string = False

                    for part in parts:
                        try:
                            float(part)
                            cell_has_numeric = True
                        except ValueError:
                            cell_has_string = True

                    if cell_has_numeric and cell_has_string:
                        cell_warnings.append(
                            f"Row {idx + 2}, Column '{col}': Cell '{cell_value}' contains mixed types "
                            f"(both numeric and non-numeric values in semicolon-separated cell). "
                            f"This column will be treated as a string column."
                        )
                        has_string = True
                    elif cell_has_numeric:
                        has_numeric = True
                    else:
                        has_string = True
                else:
                    # Single value cell
                    try:
                        float(cell_value)
                        has_numeric = True
                    except (ValueError, TypeError):
                        has_string = True

            if has_numeric and has_string:
                mixed_type_columns.append(col)

        return mixed_type_columns, cell_warnings

    def _generate_mixed_type_warning(self, mixed_columns: list[str]) -> str | None:
        """
        Generate warning about mixed-type columns.
        """
        if not mixed_columns:
            return None

        return (
            "Clarification on mixed types columns: "
            "Columns are treated as string by DivBase if they contain a mix of numeric and non-numeric values "
            "or numeric-looking values with extra characters (for example commas, hyphens, or range-like patterns such as '1-2'). "
            "A column is only numeric if all values (including each part in semicolon-separated cells) are valid numbers. "
            "Use semicolons (;) to separate multiple numeric values. "
            "Numeric query operations (ranges, inequalities) will not be applicable to string columns."
        )

    def _classify_columns(
        self, df: pd.DataFrame, mixed_type_columns: list[str]
    ) -> tuple[list[str], list[str], list[str]]:
        """
        Classify columns as numeric, string, or mixed-type.
        """
        numeric_cols = []
        string_cols = []

        for col in df.columns:
            if col == "Sample_ID":
                continue

            if col in mixed_type_columns:
                continue

            series = df[col]

            if pd.api.types.is_numeric_dtype(series) or self.is_semicolon_separated_numeric_column(series):
                numeric_cols.append(col)
            else:
                string_cols.append(col)

        return numeric_cols, string_cols, mixed_type_columns

    def _check_for_commas(self, df: pd.DataFrame) -> list[str]:
        """
        Check for comma-separated values in any cell (warns user to use semicolons).
        """
        warnings = []

        for col in df.columns:
            series = df[col]
            # Only check string columns
            if not pd.api.types.is_string_dtype(series) and not pd.api.types.is_object_dtype(series):
                continue

            for idx, cell_value in series.items():
                if isinstance(cell_value, str) and "," in cell_value:
                    warnings.append(
                        f"Row {idx + 2}, Column '{col}': Cell contains comma. "
                        "Use semicolons (;) to separate multiple values, not commas."
                    )
                    break  # Only warn once per column

        return warnings

    def _check_for_array_notation(self, df: pd.DataFrame) -> list[str]:
        """
        Check for Python/JSON-style array notation (e.g. '[1, 2, 3]') in any cell.
        Array notation is not supported by DivBase - the column will be treated as a string.
        """
        warnings = []

        for col in df.columns:
            series = df[col]
            if not pd.api.types.is_string_dtype(series) and not pd.api.types.is_object_dtype(series):
                continue

            for idx, cell_value in series.items():
                if isinstance(cell_value, str) and cell_value.startswith("[") and cell_value.endswith("]"):
                    warnings.append(
                        f"Row {idx + 2}, Column '{col}': Cell '{cell_value}' uses array notation '[...]'. "
                        "DivBase does not support Python/JSON array notation. "
                        "This column will be treated as a string. "
                        "Use semicolons (;) to separate multiple values instead (e.g., '1;2;3')."
                    )
                    break

        return warnings

    def _validate_dimensions_match(
        self, tsv_samples: set[str], project_samples: set[str]
    ) -> tuple[list[str], list[str]]:
        """
        Validate that TSV samples match project dimensions.
        """

        # TODO consider the fact that the query route also runs _check_that_dimensions_is_up_to_date_with_VCF_files_in_bucket in tasks.py before even reaching the SharedMetadataValidator...

        errors = []
        warnings = []

        missing_from_project = tsv_samples - project_samples
        if missing_from_project:
            examples = sorted(list(missing_from_project))
            errors.append(
                f"The following samples in the TSV were not found in the DivBase project's dimensions index: {examples}. "
                "DivBase requires that all samples in the TSV file must be present in the project's dimensions index to be used for queries."
            )

        missing_from_tsv = project_samples - tsv_samples
        if missing_from_tsv:
            examples = sorted(list(missing_from_tsv))
            warnings.append(
                f"The following samples in the DivBase project's dimensions index were not found in the TSV: {examples}. "
                "This is allowed for DivBase metadata TSV files, but please be aware that these samples will not be considered when making queries with this metadata file."
            )

        return errors, warnings
