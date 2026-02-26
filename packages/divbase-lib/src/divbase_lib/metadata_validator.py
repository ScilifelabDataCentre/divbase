"""
Shared validation logic for DivBase sidecar metadata TSV files.

This file contains the single source of truth for the TSV content validation logic used by both
the CLI validator (MetadataTSVValidator) on the client-side, and the SidecarQueryManager on the server-side.

Note! Logic for the queries themselves (e.g. how filtering is handled) is not shared between the two.
This file is only for validation of the contents of the TSV file, not for query processing.
"""

import ast
import csv
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any

import pandas as pd


class ValidationCategory(Enum):
    """
    Categories for validation messages to be used by the sever-side query engine in SidecarQueryManager.
    These categories are used by SidecarQueryManager to act on the validation results by Enums instead of string-matching the error and warning messages.
    """

    # Error categories
    FILE_READ = "file_read"
    HEADER = "header"
    SAMPLE_ID_COLUMN = "sample_id_column"
    SAMPLE_ID_VALUE = "sample_id_value"
    ROW_FORMAT = "row_format"
    LIST_SYNTAX = "list_syntax"
    MIXED_TYPE = "mixed_type"
    DIMENSIONS = "dimensions"
    # Warning categories
    FORMAT = "format"
    TYPE_CLASSIFICATION = "type_classification"


@dataclass
class ValidationMessage:
    """A class to store the category and message for errors and warnings detected by SharedMetadataValidator."""

    category: ValidationCategory
    message: str


@dataclass
class MetadataValidationResult:
    """
    Dataclass to hold the results of the TSV file validation. Used by the callers of SharedMetadataValidator
    (including the client-side CLI validator and the server-side query engine).
    """

    errors: list[ValidationMessage] = field(default_factory=list)
    warnings: list[ValidationMessage] = field(default_factory=list)
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
    - Sample_ID (empty, duplicates, no list values)
    - Column type (numeric vs string, list-type multi-value cells)
    - Data format (commas, whitespace, column count)
    - List syntax validation and mixed type detection


    IMPORTANT! This class never raises errors, it collects them. This is to allow the output of the class to be compatible
    with the CLI validator on the client-side and the query engine (SidecarQueryManager) on the server side. The CLI validator
    is designed to present all errors and warnings to the user in a single pass in terminal display so that they can address all of them.
    The sever-side use the errors and warnings collected by this class to raise expections on the first error it encounters in order to
    protect the query engine from malformed TSV content.
    """

    def __init__(
        self,
        file_path: Path,
        project_samples: set[str] | None = None,
        skip_dimensions_check: bool = False,
        dimensions_sample_preview_limit: int | None = 20,
    ):
        self.file_path = file_path
        self.project_samples = project_samples
        self.skip_dimensions_check = skip_dimensions_check
        self.dimensions_sample_preview_limit = dimensions_sample_preview_limit
        self.result = MetadataValidationResult()

    def load_and_validate(self) -> MetadataValidationResult:
        """
        Main entry point to the class. Load a TSV file and call helper methods to validate it.

        After loading the TSV into a pandas DataFrame, cells that look like Python
        list literals (starting with '[') are parsed with ast.literal_eval so that
        downstream validation and querying can work with native Python lists and
        their correctly-inferred element types.
        """
        try:
            with open(self.file_path, "r", newline="", encoding="utf-8") as f:
                reader = csv.reader(f, delimiter="\t")
                rows = list(reader)

            if not rows:
                self.result.errors.append(ValidationMessage(ValidationCategory.FILE_READ, "File is empty"))
                return self.result

            # Pre-pandas checks:
            first_line = "\t".join(rows[0])
            header_errors, header_warnings = self._validate_raw_header(first_line)
            self.result.errors.extend(header_errors)
            self.result.warnings.extend(header_warnings)

            if len(rows) > 1:
                row_errors, row_warnings = self._check_row_formatting(rows)
                self.result.errors.extend(row_errors)
                self.result.warnings.extend(row_warnings)

            # Initiate Pandas dataframe from TSV and check for parsing issues:
            df = pd.read_csv(self.file_path, sep="\t", skipinitialspace=True, on_bad_lines="skip")
            df.columns = df.columns.str.lstrip("#")
            self._strip_whitespace_from_cells(df)

            list_syntax_errors = self._parse_list_cells_in_dataframe(df)
            self.result.errors.extend(list_syntax_errors)
            self.result.df = df

            # Dataframe checks:
            semicolon_warnings = self._check_for_semicolons_in_plain_string_cells(df)
            self.result.warnings.extend(semicolon_warnings)

            comma_warnings = self._check_for_commas_in_plain_string_cells(df)
            self.result.warnings.extend(comma_warnings)

            sample_id_errors, sample_id_warnings = self._validate_sample_ids(df)
            self.result.errors.extend(sample_id_errors)
            self.result.warnings.extend(sample_id_warnings)

            numeric_cols, string_cols, mixed_type_columns, cell_errors, cell_warnings = self._classify_column_type(df)
            self.result.errors.extend(cell_errors)
            self.result.warnings.extend(cell_warnings)
            self.result.mixed_type_columns = mixed_type_columns
            self.result.numeric_columns = numeric_cols
            self.result.string_columns = string_cols

            mixed_type_warning = self._generate_mixed_type_warning_clarification(mixed_type_columns)
            if mixed_type_warning:
                self.result.warnings.append(mixed_type_warning)

            if not self.skip_dimensions_check and self.project_samples is not None:
                tsv_samples = set(df["Sample_ID"].tolist()) if "Sample_ID" in df.columns else set()
                dim_errors, dim_warnings = self._validate_dimensions_match(tsv_samples, self.project_samples)
                self.result.errors.extend(dim_errors)
                self.result.warnings.extend(dim_warnings)

        except Exception as e:
            self.result.errors.append(ValidationMessage(ValidationCategory.FILE_READ, f"Failed to read file: {e}"))

        return self.result

    def _strip_whitespace_from_cells(self, df: pd.DataFrame) -> None:
        """Strip leading/trailing whitespace from all string cells in the DataFrame."""
        for col in df.columns:
            df[col] = df[col].apply(lambda x: x.strip() if isinstance(x, str) else x)

    def _parse_list_cells_in_dataframe(self, df: pd.DataFrame) -> list[ValidationMessage]:
        """
        Parse all string cells in object columns that look like Python list literals, and collect errors for cells that fail to parse.

        Only columns with dtype "object" are considered since numeric columns inferred by pandas cannot contain list strings.
        Cells that start with '[' are parsed via ast.literal_eval. On success, the cell is replaced with the parsed Python list; on failure the cell is left as-is
        and an error message is collected.

        ast.literal_eval is whitespace-insensitive within list notation:
        [3,2], [3, 2], and [ 3 , 2 ] all parse identically to [3, 2].
        """
        errors: list[ValidationMessage] = []

        for col in df.select_dtypes(include=["object"]).columns:
            for idx, cell_value in df[col].items():
                if not isinstance(cell_value, str):
                    continue
                stripped = cell_value.strip()
                if not stripped.startswith("["):
                    continue
                try:
                    parsed = ast.literal_eval(stripped)
                    if isinstance(parsed, list):
                        df.at[idx, col] = parsed
                    else:
                        errors.append(
                            ValidationMessage(
                                ValidationCategory.LIST_SYNTAX,
                                f"Row {idx + 2}, Column '{col}': Cell '{cell_value}' starts with '[' "
                                f"but parsed as {type(parsed).__name__}, not a list.",
                            )
                        )
                except (ValueError, SyntaxError):
                    errors.append(
                        ValidationMessage(
                            ValidationCategory.LIST_SYNTAX,
                            f"Row {idx + 2}, Column '{col}': Cell '{cell_value}' has invalid "
                            "Python list syntax. Multi-value cells must use valid Python list "
                            'notation, e.g. [1, 2, 3] or ["a", "b"].',
                        )
                    )

        return errors

    def _validate_raw_header(self, header_line: str) -> tuple[list[ValidationMessage], list[ValidationMessage]]:
        """
        Validate the raw header line before pandas processing.

        Attempts to catch issues that pandas would silently fix (like duplicate columns).
        """
        errors: list[ValidationMessage] = []

        raw_columns = header_line.split("\t")
        cleaned_columns = [col.lstrip("#") for col in raw_columns]

        empty_columns = [i + 1 for i, col in enumerate(cleaned_columns) if not col.strip()]
        if empty_columns:
            errors.append(
                ValidationMessage(
                    ValidationCategory.HEADER,
                    f"Empty column name(s) found at position(s): {empty_columns}. All columns must have a non-empty name.",
                )
            )

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
                ValidationMessage(
                    ValidationCategory.HEADER,
                    f"Duplicate column names found: {duplicate_columns}. "
                    "Each column name must be unique in the metadata file.",
                )
            )

        if raw_columns and raw_columns[0] != "#Sample_ID":
            errors.append(
                ValidationMessage(
                    ValidationCategory.SAMPLE_ID_COLUMN,
                    f"First column must be named '#Sample_ID', found: '{raw_columns[0]}'",
                )
            )

        return errors, []

    def _check_row_formatting(self, rows: list[list[str]]) -> tuple[list[ValidationMessage], list[ValidationMessage]]:
        """
        Check for row-level formatting issues that pandas might handle silently.
        """
        errors: list[ValidationMessage] = []
        warnings: list[ValidationMessage] = []

        header = rows[0]
        num_columns = len(header)

        for row_num, row in enumerate(rows[1:], start=2):  # Start at row 2 (i.e. skip header)
            if len(row) != num_columns:
                sample_hint = f" (Sample_ID: '{row[0]}')" if row else ""
                errors.append(
                    ValidationMessage(
                        ValidationCategory.ROW_FORMAT,
                        f"Row {row_num}: Expected {num_columns} tab-separated columns from reading the header, "
                        f"found {len(row)}{sample_hint}. "
                        "Check that all cells in the TSV are separated by tabs (not spaces).",
                    )
                )
                continue

            for col_idx, cell in enumerate(row):
                if cell != cell.strip():
                    col_name = header[col_idx]
                    warnings.append(
                        ValidationMessage(
                            ValidationCategory.FORMAT,
                            f"Row {row_num}, Column '{col_name}': Cell has leading or trailing whitespace "
                            "(this is allowed, but note that they will be stripped by DivBase server when the TSV is used for queries)",
                        )
                    )

                if col_idx == 0 and not cell.strip():
                    errors.append(
                        ValidationMessage(
                            ValidationCategory.SAMPLE_ID_VALUE,
                            f"Row {row_num}: Sample_ID is empty",
                        )
                    )

        return errors, warnings

    def _validate_sample_ids(self, df: pd.DataFrame) -> tuple[list[ValidationMessage], list[ValidationMessage]]:
        """
        Validate Sample_ID column in the DataFrame.

        The Sample_ID column must be present in the TSV, it has to be the first column, and it must contain non-empty, unique values.
        Further more, it must have a single value per row. List values (Python list notation like ["S1", "S2"]) are not allowed in Sample_ID.
        """
        errors: list[ValidationMessage] = []

        if "Sample_ID" not in df.columns:
            errors.append(
                ValidationMessage(
                    ValidationCategory.SAMPLE_ID_COLUMN,
                    "The 'Sample_ID' column is required in the metadata file.",
                )
            )
            return errors, []

        if df["Sample_ID"].isna().any() or (df["Sample_ID"] == "").any():
            errors.append(
                ValidationMessage(
                    ValidationCategory.SAMPLE_ID_VALUE,
                    "Sample_ID column contains empty or missing values. All rows must have a valid Sample_ID.",
                )
            )

        if df["Sample_ID"].duplicated().any():
            duplicates = df[df["Sample_ID"].duplicated()]["Sample_ID"].tolist()
            errors.append(
                ValidationMessage(
                    ValidationCategory.SAMPLE_ID_VALUE,
                    f"Duplicate Sample_IDs found: {duplicates}. Each Sample_ID must be unique.",
                )
            )

        list_sample_ids = [
            sid
            for sid in df["Sample_ID"].dropna()
            if isinstance(sid, list)
            or (isinstance(sid, str) and sid.strip().startswith("[") and sid.strip().endswith("]"))
        ]
        if list_sample_ids:
            errors.append(
                ValidationMessage(
                    ValidationCategory.SAMPLE_ID_VALUE,
                    f"Sample_ID column contains list values: {list_sample_ids}. "
                    "Sample_ID must contain only one value per row (list notation is not allowed).",
                )
            )

        return errors, []

    def parse_cell_value(self, cell_value) -> Any:
        """
        Parse a single cell value. If the string representation starts with '[', it must be a valid Python list literal (parsed via ast.literal_eval).
        Non-list cells are returned as-is.

        Raises ValueError if a cell looks like a list (starts with '[') but cannot be parsed by ast.literal_eval.
        """
        if pd.isna(cell_value):
            return cell_value

        cell_str = str(cell_value).strip()
        if cell_str.startswith("["):
            try:
                parsed = ast.literal_eval(cell_str)
            except (ValueError, SyntaxError) as exc:
                raise ValueError(
                    f"Invalid Python list syntax: '{cell_str}'. "
                    "Multi-value cells must use valid Python list notation, "
                    'e.g. [1, 2, 3] or ["a", "b"].'
                ) from exc
            if not isinstance(parsed, list):
                raise ValueError(
                    f"Cell '{cell_str}' parsed successfully but is not a list (got {type(parsed).__name__}). "
                    "Multi-value cells must be Python lists."
                )
            return parsed

        return cell_value

    def _classify_column_type(
        self, df: pd.DataFrame
    ) -> tuple[list[str], list[str], list[str], list[ValidationMessage], list[ValidationMessage]]:
        """
        Classify every user-defined column as numeric, string, or mixed-type in a single pass over the data.

        For each column, every non-null cell is examined once:
        - Single-value cells: numeric if ``float()`` succeeds, otherwise string.
        - Multi-value cells (Python lists): numeric if all elements are int/float, string if all are strings.
          A list with mixed element types (e.g. [1, "two"]) is a hard error because ast.literal_eval preserves the exact types
          the user wrote, so mixed types indicate an explicit mistake.

        After scanning all cells in a column, classify the column type:
        - All cells are numeric -> numeric column
        - All cells are string  -> string column
        - Contains both numeric and string cells -> mixed-type column: treat as string and send per-cell warnings to communicate the ambiguities to the user.
        - All cells are Null -> numeric (to match Pandas default for NaN-only columns)
        """
        numeric_cols: list[str] = []
        string_cols: list[str] = []
        mixed_type_columns: list[str] = []
        cell_errors: list[ValidationMessage] = []
        cell_warnings: list[ValidationMessage] = []

        for col in df.columns:
            if col == "Sample_ID":
                continue

            series = df[col]
            non_null_values = series.dropna()

            if len(non_null_values) == 0:
                numeric_cols.append(col)
                continue

            has_numeric = False
            has_string = False
            numeric_cells: list[tuple[int, Any]] = []
            string_cells: list[tuple[int, Any]] = []

            for idx, cell_value in non_null_values.items():
                if isinstance(cell_value, list):
                    cell_has_numeric = False
                    cell_has_string = False

                    for element in cell_value:
                        if isinstance(element, (int, float)):
                            cell_has_numeric = True
                        else:
                            cell_has_string = True

                    if cell_has_numeric and cell_has_string:
                        cell_errors.append(
                            ValidationMessage(
                                ValidationCategory.MIXED_TYPE,
                                f"Row {idx + 2}, Column '{col}': List cell {cell_value} contains "
                                f"mixed element types (both numeric and string values). "
                                f"All elements in a list must be the same type. Use "
                                f"either all numbers (e.g. [1, 2, 3]) or all strings "
                                f'(e.g. ["a", "b", "c"]).',
                            )
                        )
                        has_string = True
                        string_cells.append((idx, cell_value))
                    elif cell_has_numeric:
                        has_numeric = True
                        numeric_cells.append((idx, cell_value))
                    else:
                        has_string = True
                        string_cells.append((idx, cell_value))
                else:
                    try:
                        float(cell_value)
                        has_numeric = True
                        numeric_cells.append((idx, cell_value))
                    except (ValueError, TypeError):
                        has_string = True
                        string_cells.append((idx, cell_value))

            if has_numeric and has_string:
                mixed_type_columns.append(col)
                minority = string_cells if len(string_cells) <= len(numeric_cells) else numeric_cells
                minority_type = "non-numeric" if minority is string_cells else "numeric"
                for idx, val in minority:
                    cell_warnings.append(
                        ValidationMessage(
                            ValidationCategory.TYPE_CLASSIFICATION,
                            f"Row {idx + 2}, Column '{col}': Cell '{val}' is {minority_type} "
                            f"in a column that contains both numeric and non-numeric values. "
                            f"This column will be treated as string type.",
                        )
                    )
            elif has_numeric:
                numeric_cols.append(col)
            else:
                string_cols.append(col)

        return numeric_cols, string_cols, mixed_type_columns, cell_errors, cell_warnings

    def _generate_mixed_type_warning_clarification(self, mixed_columns: list[str]) -> ValidationMessage | None:
        """
        Generate clarification warning about mixed-type columns. This tells users that mixed-type columns are treated as string in Divbase.
        """
        if not mixed_columns:
            return None

        return ValidationMessage(
            ValidationCategory.TYPE_CLASSIFICATION,
            "Clarification on mixed types columns: "
            "Columns are treated as string by DivBase if they contain a mix of numeric and non-numeric values "
            'or list cells with mixed element types (for example [1, "two"]). '
            "A column is only numeric if all values (including each element in list cells) are valid numbers. "
            "Use Python list notation (e.g. [1, 2, 3]) for multi-value cells. "
            "Numeric query operations (ranges, inequalities) will not be applicable to string columns.",
        )

    def _check_for_semicolons_in_plain_string_cells(self, df: pd.DataFrame) -> list[ValidationMessage]:
        """
        Check for semicolons in plain string (single-value) cells. This is to reduce user confusion with the DivBase query filter syntax,
        which uses semicolons to separate key:value pairs. For example: divbase query tsv "Area:North;Population:1"
        filters for rows where Area is "North" AND Population is 1.

        A TSV cell containing a semicolon (e.g. "2;4") will be treated as a plain string value and cannot be matched via the query syntax since
        the query parser will split on the semicolon. If the user intended multiple values, they should use Python list notation instead.
        """
        warnings: list[ValidationMessage] = []

        for col in df.columns:
            if col == "Sample_ID":
                continue
            series = df[col]
            if not pd.api.types.is_string_dtype(series) and not pd.api.types.is_object_dtype(series):
                continue

            for idx, cell_value in series.items():
                if not isinstance(cell_value, str):
                    continue
                if ";" in cell_value:
                    warnings.append(
                        ValidationMessage(
                            ValidationCategory.FORMAT,
                            f"Row {idx + 2}, Column '{col}': Cell '{cell_value}' contains a semicolon. "
                            "If you intended multiple values, use Python list notation instead "
                            '(e.g. [1, 2] or ["a", "b"]). '
                            "If the semicolon is intentional, note that DivBase query syntax uses "
                            "semicolons to separate filter key:value pairs, so this exact cell value "
                            "cannot be matched via queries.",
                        )
                    )
                    break

        return warnings

    def _check_for_commas_in_plain_string_cells(self, df: pd.DataFrame) -> list[ValidationMessage]:
        """
        Check for commas in plain string cells (in single-value cells, not in multi-values list cells) and warn users about the
        ambiguity that might cause for DivBase filtering since the metadata query filter syntax uses commas to separate filter values.
        For example: divbase query tsv "Area:North,South" filters for rows where Area is "North" OR "South".

        A TSV cell containing the literal string "North,South" would not match that query, because the query parser splits on commas.
        This helper method warns users about this.

        Commas inside list notation (e.g. ["North", "South"]) are fine since they are parsed as lists.
        """
        warnings: list[ValidationMessage] = []

        for col in df.columns:
            series = df[col]
            if not pd.api.types.is_string_dtype(series) and not pd.api.types.is_object_dtype(series):
                continue

            for idx, cell_value in series.items():
                if not isinstance(cell_value, str):
                    continue
                if "," in cell_value:
                    warnings.append(
                        ValidationMessage(
                            ValidationCategory.FORMAT,
                            f"Row {idx + 2}, Column '{col}': Cell '{cell_value}' contains a comma. "
                            "If you intended multiple values, use Python list notation instead "
                            '(e.g. [1, 2] or ["a", "b"]). '
                            "If the comma is intentional, note that DivBase query syntax uses "
                            "commas to separate filter values. To query for this exact string, "
                            "enclose the value in double quotes in your filter (e.g. "
                            'Area:"North,South").',
                        )
                    )
                    break

        return warnings

    def _validate_dimensions_match(
        self, tsv_samples: set[str], project_samples: set[str]
    ) -> tuple[list[ValidationMessage], list[ValidationMessage]]:
        """
        Validate that TSV samples match project dimensions.
        """

        # TODO consider the fact that the query route also runs _check_that_dimensions_is_up_to_date_with_VCF_files_in_bucket in tasks.py before even reaching the SharedMetadataValidator...

        errors: list[ValidationMessage] = []
        warnings: list[ValidationMessage] = []

        missing_from_project = tsv_samples - project_samples
        if missing_from_project:
            examples, was_truncated = self._format_sample_name_preview(
                sample_names=missing_from_project, preview_limit=self.dimensions_sample_preview_limit
            )
            full_output_hint = f" {self._build_full_sample_output_hint()}" if was_truncated else ""
            errors.append(
                ValidationMessage(
                    ValidationCategory.DIMENSIONS,
                    f"The following samples in the TSV were not found in the DivBase project's dimensions index ({examples}). "
                    "DivBase requires that all samples in the TSV file must be present in the project's dimensions index to be used for queries."
                    f"{full_output_hint}",
                )
            )

        missing_from_tsv = project_samples - tsv_samples
        if missing_from_tsv:
            examples, was_truncated = self._format_sample_name_preview(
                sample_names=missing_from_tsv, preview_limit=self.dimensions_sample_preview_limit
            )
            full_output_hint = f" {self._build_full_sample_output_hint()}" if was_truncated else ""
            warnings.append(
                ValidationMessage(
                    ValidationCategory.DIMENSIONS,
                    f"The following samples in the DivBase project's dimensions index were not found in the TSV ({examples}). "
                    "This is allowed for DivBase metadata TSV files, but please be aware that these samples will not be considered when making queries with this metadata file."
                    f"{full_output_hint}",
                )
            )

        return errors, warnings

    def _format_sample_name_preview(self, sample_names: set[str], preview_limit: int | None = 20) -> tuple[str, bool]:
        """
        Build a compact sample-name summary for terminal messages. Truncates list of sample names if it exceeds the preview limit (20 by default).

        Returns only a full list when the list is small, otherwise includes a count and
        a preview of the first 20 samples (default value) to avoid overwhelming CLI output.
        """
        sorted_names = sorted(sample_names)
        total = len(sorted_names)
        if preview_limit is None:
            return f"count: {total}, samples: {sorted_names}", False

        if total <= preview_limit:
            return f"count: {total}, samples: {sorted_names}", False

        preview = sorted_names[:preview_limit]
        return f"count: {total}, showing first {preview_limit}: {preview}", True

    def _build_full_sample_output_hint(self) -> str:
        """
        When the sample mismatch list is truncated (>20 for default settings), give users a hint on how they can see the full sample-name mismatch list.
        """
        return (
            "To view the full list of mismatched samples, run: "
            "divbase-cli dimensions validate-metadata-file <TSV_FILE_NAME> --project <PROJECT_NAME> "
            "--full-sample-mismatch-names."
        )
