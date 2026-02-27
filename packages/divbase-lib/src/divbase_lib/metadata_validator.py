"""
Shared validation logic for DivBase sidecar metadata TSV files.

This file contains the single source of truth for the TSV content validation logic used by both
the CLI validator (ClientSideMetadataTSVValidator) on the client-side, and the SidecarQueryManager on the server-side.

Note! Logic for the queries themselves (e.g. how filtering is handled) is not shared between the two.
This file is only for validation of the contents of the TSV file, not for query processing.
"""

import ast
import csv
from collections import defaultdict
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
class ValidationStats:
    """A class to store any relevant statistics collected during validation. Used on the CLI side."""

    total_columns: int = 0
    user_defined_columns: int = 0
    samples_in_tsv: int = 0
    samples_matching_project: int = 0
    total_project_samples: int = 0
    numeric_column_count: int = 0
    string_column_count: int = 0
    mixed_type_column_count: int = 0
    empty_cells_per_column: dict[str, int] = field(default_factory=dict)
    has_multi_values: bool = False


@dataclass
class MetadataValidationResult:
    """
    Dataclass to hold the results of the TSV file validation. Used by the callers of SharedMetadataValidator
    (including the client-side CLI validator and the server-side query engine).
    """

    errors: list[ValidationMessage] = field(default_factory=list)
    warnings: list[ValidationMessage] = field(default_factory=list)
    stats: ValidationStats = field(default_factory=ValidationStats)
    df: pd.DataFrame | None = None
    mixed_type_columns: list[str] = field(default_factory=list)
    numeric_columns: list[str] = field(default_factory=list)
    string_columns: list[str] = field(default_factory=list)


class SharedMetadataValidator:
    """
    Core validation logic for DivBase sidecar metadata TSV files. Shared between client-side ClientSideMetadataTSVValidator
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
        self.tsv_samples: set[str] = set()
        self.skip_dimensions_check = skip_dimensions_check
        self.dimensions_sample_preview_limit = dimensions_sample_preview_limit
        self._cells_with_hard_errors: set[tuple[int, str]] = set()
        self.result = MetadataValidationResult()

    def load_and_validate(self) -> MetadataValidationResult:
        """
        Main entry point to the class. Load a TSV file and call helper methods to validate it.

        After loading the TSV into a pandas DataFrame, cells that look like Python
        list literals (starting with '[') are parsed with ast.literal_eval so that
        downstream validation and querying can work with native Python lists and
        their correctly-inferred element types.
        """
        # Pre-pandas checks:
        try:
            with open(self.file_path, "r", newline="", encoding="utf-8") as f:
                reader = csv.reader(f, delimiter="\t")
                rows = list(reader)
        except Exception as e:
            self.result.errors.append(
                ValidationMessage(ValidationCategory.FILE_READ, f"Failed to read file (TSV read stage): {e}")
            )
            return self.result

        if not rows:
            self.result.errors.append(ValidationMessage(ValidationCategory.FILE_READ, "File is empty"))
            return self.result

        first_line = "\t".join(rows[0])
        header_errors, header_warnings = self._validate_raw_header(first_line)
        self.result.errors.extend(header_errors)
        self.result.warnings.extend(header_warnings)

        if len(rows) > 1:
            row_errors, row_warnings = self._check_row_formatting(rows)
            self.result.errors.extend(row_errors)
            self.result.warnings.extend(row_warnings)

        # Initiate Pandas dataframe from TSV and check for parsing issues:
        try:
            df = pd.read_csv(self.file_path, sep="\t", skipinitialspace=True, on_bad_lines="skip")
            df.columns = df.columns.str.lstrip("#")
            self._strip_whitespace_from_cells(df)
        except Exception as e:
            self.result.errors.append(
                ValidationMessage(ValidationCategory.FILE_READ, f"Failed to read file (pandas parse stage): {e}")
            )
            return self.result

        # Dataframe checks:
        try:
            list_syntax_errors, cells_with_hard_errors = self._parse_list_cells_in_dataframe(df)
            self.result.errors.extend(list_syntax_errors)
            self._cells_with_hard_errors = cells_with_hard_errors
            self.result.df = df

            semicolon_warnings = self._check_for_semicolons_in_plain_string_cells(df)
            self.result.warnings.extend(semicolon_warnings)

            comma_warnings = self._check_for_commas_in_plain_string_cells(df)
            self.result.warnings.extend(comma_warnings)

            sample_id_errors, sample_id_warnings = self._validate_sample_ids(df)
            self.result.errors.extend(sample_id_errors)
            self.result.warnings.extend(sample_id_warnings)

            numeric_cols, string_cols, mixed_type_cols, cell_errors, cell_warnings = self._classify_column_type(df)
            self.result.errors.extend(cell_errors)
            self.result.warnings.extend(cell_warnings)
            self.result.mixed_type_columns = mixed_type_cols
            self.result.numeric_columns = numeric_cols
            self.result.string_columns = string_cols

            self.tsv_samples = self._extract_sample_ids_for_dimensions(df)

            if not self.skip_dimensions_check and self.project_samples is not None:
                dim_errors, dim_warnings = self._validate_dimensions_match(self.tsv_samples, self.project_samples)
                self.result.errors.extend(dim_errors)
                self.result.warnings.extend(dim_warnings)
        except Exception as e:
            self.result.errors.append(ValidationMessage(ValidationCategory.FILE_READ, f"Validation failed: {e}"))

        self.result.stats = self._collect_statistics(df)

        return self.result

    def _strip_whitespace_from_cells(self, df: pd.DataFrame) -> None:
        """Strip leading/trailing whitespace from all string cells in the DataFrame."""
        for col in df.columns:
            df[col] = df[col].apply(lambda x: x.strip() if isinstance(x, str) else x)

    def _parse_list_cells_in_dataframe(self, df: pd.DataFrame) -> tuple[list[ValidationMessage], set[tuple[int, str]]]:
        """
        Parse all string cells in object columns that look like Python list literals, and collect errors for cells that fail to parse.

        Only columns with dtype "object" are considered since numeric columns inferred by pandas cannot contain list strings.
        Cells that start with '[' are parsed via ast.literal_eval. On success, the cell is replaced with the parsed Python list; on failure the cell is left as-is
        and an error message is collected.

        ast.literal_eval is whitespace-insensitive within list notation:
        [3,2], [3, 2], and [ 3 , 2 ] all parse identically to [3, 2].
        """
        errors: list[ValidationMessage] = []
        cells_with_hard_errors: set[tuple[int, str]] = set()
        mixed_type_cells_by_column: defaultdict[str, list[tuple[int, Any]]] = defaultdict(
            list
        )  # Automatically initializes empty lists for new columns
        list_syntax_errors_by_column: defaultdict[str, list[tuple[int, Any]]] = defaultdict(list)

        # Detect errors on a cell basis
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
                        # Also reject lists that contain numeric and non-numeric types (e.g. [1, "two"]).
                        has_numeric = any(isinstance(element, (int, float)) for element in parsed)
                        has_non_numeric = any(not isinstance(element, (int, float)) for element in parsed)
                        if has_numeric and has_non_numeric:
                            mixed_type_cells_by_column[col].append((idx + 2, parsed))
                            cells_with_hard_errors.add((idx, col))
                            continue
                        df.at[idx, col] = parsed
                    else:
                        list_syntax_errors_by_column[col].append((idx + 2, cell_value))
                        cells_with_hard_errors.add((idx, col))
                except (ValueError, SyntaxError):
                    list_syntax_errors_by_column[col].append((idx + 2, cell_value))
                    cells_with_hard_errors.add((idx, col))

        # Collect cell value errors by column for brevity
        for col, mixed_cells in mixed_type_cells_by_column.items():
            preview_text, was_truncated = self._format_row_value_preview(
                values=mixed_cells,
                preview_limit=self.dimensions_sample_preview_limit,
            )
            full_output_hint = f" {self._build_full_sample_output_hint()}" if was_truncated else ""
            errors.append(
                ValidationMessage(
                    ValidationCategory.MIXED_TYPE,
                    f"Column '{col}': Found {len(mixed_cells)} cell(s) with mixed element types (both numeric and string values) in lists. "
                    'All elements in a list must be the same type. Use either all numbers (e.g. [1, 2, 3]) or all strings (e.g. ["a", "b", "c"]). '
                    f"Affected cells: {preview_text}.{full_output_hint}",
                )
            )
        for col, error_cells in list_syntax_errors_by_column.items():
            preview_text, was_truncated = self._format_row_value_preview(
                values=error_cells,
                preview_limit=self.dimensions_sample_preview_limit,
            )
            full_output_hint = f" {self._build_full_sample_output_hint()}" if was_truncated else ""
            errors.append(
                ValidationMessage(
                    ValidationCategory.LIST_SYNTAX,
                    f"Column '{col}': Found {len(error_cells)} cell(s) with invalid list syntax or not parsed as list. "
                    'Multi-value cells must use valid Python list notation, e.g. [1, 2, 3] or ["a", "b"]. '
                    f"Affected cells: {preview_text}.{full_output_hint}",
                )
            )

        return errors, cells_with_hard_errors

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
        whitespace_cells_by_column: defaultdict[str, list[tuple[int, str]]] = defaultdict(list)

        header = rows[0]
        num_columns = len(header)

        for row_num, row in enumerate(rows[1:], start=2):  # Start at row 2 (i.e. skip header)
            # Report tab separation issues by row
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
                    whitespace_cells_by_column[col_name].append((row_num, cell))

        for col_name, whitespace_cells in whitespace_cells_by_column.items():
            preview_text, was_truncated = self._format_row_value_preview(
                values=whitespace_cells, preview_limit=self.dimensions_sample_preview_limit
            )
            full_output_hint = f" {self._build_full_sample_output_hint()}" if was_truncated else ""
            warnings.append(
                ValidationMessage(
                    ValidationCategory.FORMAT,
                    f"Column '{col_name}': Found {len(whitespace_cells)} cell(s) with leading or trailing whitespace "
                    "(this is allowed, but note that they will be stripped by DivBase server when the TSV is used for queries). "
                    f"Affected cells: {preview_text}.{full_output_hint}",
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
        empty_sample_id_rows = []

        if "Sample_ID" not in df.columns:
            errors.append(
                ValidationMessage(
                    ValidationCategory.SAMPLE_ID_COLUMN,
                    "The 'Sample_ID' column is required in the metadata file.",
                )
            )
            return errors, []

        for idx, val in df["Sample_ID"].items():
            # Avoid pd.isna(list/array) ambiguity and treat only missing/empty scalar values as missing Sample_ID.
            if val is None or (isinstance(val, str) and val == "") or (not isinstance(val, list) and pd.isna(val)):
                empty_sample_id_rows.append(
                    (idx + 2, "<empty>")
                )  # If val is used here, it will display pandas 'nan' to the user. <empty> is a clearer message
        if empty_sample_id_rows:
            preview_text, was_truncated = self._format_row_value_preview(
                values=empty_sample_id_rows,
                preview_limit=self.dimensions_sample_preview_limit,
            )
            full_output_hint = f" {self._build_full_sample_output_hint()}" if was_truncated else ""
            errors.append(
                ValidationMessage(
                    ValidationCategory.SAMPLE_ID_VALUE,
                    f"Sample_ID is empty or missing in {len(empty_sample_id_rows)} row(s). "
                    f"Affected rows: {preview_text}.{full_output_hint} All rows must have a valid Sample_ID.",
                )
            )

        duplicate_rows_by_value: defaultdict[str, list[tuple[int, str]]] = defaultdict(list)
        for idx, val in df["Sample_ID"].items():
            if isinstance(val, list):
                continue
            if val is None or (isinstance(val, str) and val == "") or (not isinstance(val, list) and pd.isna(val)):
                continue
            sample_id_str = str(val)
            duplicate_rows_by_value[sample_id_str].append((idx + 2, sample_id_str))

        for sample_id, rows in duplicate_rows_by_value.items():
            if len(rows) <= 1:
                continue
            preview_text, was_truncated = self._format_row_value_preview(
                values=rows,
                preview_limit=self.dimensions_sample_preview_limit,
            )
            full_output_hint = f" {self._build_full_sample_output_hint()}" if was_truncated else ""
            errors.append(
                ValidationMessage(
                    ValidationCategory.SAMPLE_ID_VALUE,
                    f"Duplicate Sample_IDs found: '{sample_id}' appears in {len(rows)} row(s). "
                    f"Affected rows: {preview_text}.{full_output_hint} Each Sample_ID must be unique.",
                )
            )

        list_syntax_in_sample_ids = [
            sid
            for sid in df["Sample_ID"].dropna()
            if isinstance(sid, list)
            or (isinstance(sid, str) and sid.strip().startswith("[") and sid.strip().endswith("]"))
        ]
        if list_syntax_in_sample_ids:
            errors.append(
                ValidationMessage(
                    ValidationCategory.SAMPLE_ID_VALUE,
                    f"Sample_ID column contains list values: {list_syntax_in_sample_ids}. "
                    "Sample_ID must contain only one value per row (list notation is not allowed).",
                )
            )

        return errors, []

    def _extract_sample_ids_for_dimensions(self, df: pd.DataFrame) -> set[str]:
        """
        Extract sample IDs for dimensions matching in a way that avoid pandas type-related exceptions that can be cryptic to users.

        Skips missing values and list-valued Sample_ID cells (e.g. ['test_duplicate', 'test_duplicate2']) which are handled by
        Sample_ID validation errors elsewhere.
        """
        if "Sample_ID" not in df.columns:
            return set()

        sample_ids: set[str] = set()
        for val in df["Sample_ID"].tolist():
            if isinstance(val, list):
                continue
            if val is None or (isinstance(val, str) and val == "") or (not isinstance(val, list) and pd.isna(val)):
                continue
            sample_ids.add(str(val))
        return sample_ids

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
        - Contains both numeric and string cells -> mixed-type column: treat as string and send a single warning per column summarizing all type issues.
        - All cells are Null -> numeric (to match Pandas default for NaN-only columns)
        """
        numeric_cols: list[str] = []
        string_cols: list[str] = []
        mixed_type_cols: list[str] = []
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
            string_cells: list[tuple[int, Any]] = []  # true string cells only
            mixed_type_cells: list[tuple[int, Any]] = []

            for idx, cell_value in non_null_values.items():
                if (idx, col) in self._cells_with_hard_errors:
                    # Do not classify cells that already have hard validation errors.
                    continue

                if isinstance(cell_value, list):
                    cell_has_numeric = False
                    cell_has_string = False

                    for element in cell_value:
                        if isinstance(element, (int, float)):
                            cell_has_numeric = True
                        else:
                            cell_has_string = True

                    if cell_has_numeric and cell_has_string:
                        mixed_type_cells.append((idx + 2, cell_value))
                        has_string = True
                        # Do NOT add to string_cells; mixed-type lists are reported separately
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

            # Report mixed-type warning per column. There are two cases: cells with different types, and multi-value cells that contain mixed types in the list.
            if (has_numeric and has_string) or mixed_type_cells:
                mixed_type_cols.append(col)

                n_numeric = len(numeric_cells)
                n_string = len(string_cells)
                n_mixed = len(mixed_type_cells)

                mixed_summary = f", {n_mixed} mixed-type list cells" if n_mixed > 0 else ""

                summary = f"Column '{col}': This column contains mixed-type cells ({n_numeric} numeric, {n_string} string type cells{mixed_summary}), and will be treated as string type. Numeric query operations (ranges, inequalities) are not applicable to string columns."
                message_parts = [summary]
                hard_error_cells_in_column = sum(1 for _, error_col in self._cells_with_hard_errors if error_col == col)
                if hard_error_cells_in_column > 0:
                    message_parts.append(
                        f"{hard_error_cells_in_column} cell(s) in this column have hard list-format/type errors and are reported separately under ERRORS."
                    )

                # Only return the cell type that is in minority in the column. If equal, don't print this message
                if n_numeric > 0 and n_string > 0 and n_numeric != n_string:
                    minority_cells = string_cells if n_string <= n_numeric else numeric_cells
                    minority_type = "string" if minority_cells is string_cells else "numeric"
                    preview_text, was_truncated = self._format_row_value_preview(
                        values=[(idx + 2, val) for idx, val in minority_cells],
                        preview_limit=self.dimensions_sample_preview_limit,
                    )
                    full_output_hint = f" {self._build_full_sample_output_hint()}" if was_truncated else ""
                    message_parts.append(
                        f"The {minority_type} cells are in minority and are: {preview_text}.{full_output_hint}"
                    )

                # Always return cell values with mixed-type list cells
                if n_mixed > 0:
                    preview_text, was_truncated = self._format_row_value_preview(
                        values=mixed_type_cells,
                        preview_limit=self.dimensions_sample_preview_limit,
                    )
                    full_output_hint = f" {self._build_full_sample_output_hint()}" if was_truncated else ""
                    message_parts.append(f"Mixed-type list cells: {preview_text}.{full_output_hint}")

                cell_warnings.append(
                    ValidationMessage(
                        ValidationCategory.TYPE_CLASSIFICATION,
                        " ".join(message_parts),
                    )
                )
            elif has_numeric:
                numeric_cols.append(col)
            else:
                string_cols.append(col)

        return numeric_cols, string_cols, mixed_type_cols, cell_errors, cell_warnings

    def _check_for_semicolons_in_plain_string_cells(self, df: pd.DataFrame) -> list[ValidationMessage]:
        """
        Check for semicolons in plain string (single-value) cells. This is to reduce user confusion with the DivBase query filter syntax,
        which uses semicolons to separate key:value pairs. For example: divbase query tsv "Area:North;Population:1"
        filters for rows where Area is "North" AND Population is 1.

        A TSV cell containing a semicolon (e.g. "2;4") will be treated as a plain string value and cannot be matched via the query syntax since
        the query parser will split on the semicolon. If the user intended multiple values, they should use Python list notation instead.
        """
        warnings: list[ValidationMessage] = []
        semicolon_cells_by_column: defaultdict[str, list[tuple[int, str]]] = defaultdict(list)

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
                    semicolon_cells_by_column[col].append((idx + 2, cell_value))
                    break

        for col, warning_cells in semicolon_cells_by_column.items():
            preview_text, was_truncated = self._format_row_value_preview(
                values=warning_cells,
                preview_limit=self.dimensions_sample_preview_limit,
            )
            full_output_hint = f" {self._build_full_sample_output_hint()}" if was_truncated else ""
            warnings.append(
                ValidationMessage(
                    ValidationCategory.FORMAT,
                    f"Column '{col}': Found {len(warning_cells)} cell(s) containing a semicolon. "
                    "Note that DivBase query syntax uses semicolons to separate filter key:value pairs, "
                    "so this exact cell value cannot be matched via queries. "
                    f"Affected cells: {preview_text}.{full_output_hint}",
                )
            )

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
        comma_cells_by_column: defaultdict[str, list[tuple[int, str]]] = defaultdict(list)

        for col in df.columns:
            series = df[col]
            if not pd.api.types.is_string_dtype(series) and not pd.api.types.is_object_dtype(series):
                continue

            for idx, cell_value in series.items():
                if not isinstance(cell_value, str):
                    continue
                if "," in cell_value:
                    comma_cells_by_column[col].append((idx + 2, cell_value))

        for col, warning_cells in comma_cells_by_column.items():
            preview_text, was_truncated = self._format_row_value_preview(
                values=warning_cells,
                preview_limit=self.dimensions_sample_preview_limit,
            )
            full_output_hint = f" {self._build_full_sample_output_hint()}" if was_truncated else ""
            warnings.append(
                ValidationMessage(
                    ValidationCategory.FORMAT,
                    f"Column '{col}': Found {len(warning_cells)} cell(s) containing a comma. "
                    "If you intended to describe multiple values, use Python list notation instead "
                    '(e.g. [1, 2] or ["a", "b"]). '
                    "If the comma is intentional, note that DivBase query syntax uses "
                    "commas to separate filter values. To query for this exact string, "
                    "enclose the value in double quotes in your filter (e.g. "
                    'Area:"North,South").'
                    f"Affected cells: {preview_text}.{full_output_hint}",
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
        sorted_names = sorted(sample_names, key=lambda sample_name: str(sample_name))
        total = len(sorted_names)
        if preview_limit is None:
            return f"count: {total}, samples: {sorted_names}", False

        if total <= preview_limit:
            return f"count: {total}, samples: {sorted_names}", False

        preview = sorted_names[:preview_limit]
        return f"count: {total}, showing first {preview_limit}: {preview}", True

    def _format_row_value_preview(
        self,
        values: list[tuple[int, Any]],
        preview_limit: int | None = 20,
    ) -> tuple[str, bool]:
        """
        Collect, and when necessary, truncate, row numbers and cell values for error/warning messages about specific cells.
        """
        if preview_limit is None:
            return ", ".join(f"Row {row_num} ('{value}')" for row_num, value in values), False

        if len(values) <= preview_limit:
            return ", ".join(f"Row {row_num} ('{value}')" for row_num, value in values), False

        preview = ", ".join(f"Row {row_num} ('{value}')" for row_num, value in values[:preview_limit])
        return f"{preview}, ... ({len(values) - preview_limit} more)", True

    def _build_full_sample_output_hint(self) -> str:
        """
        When a validator list is truncated in the error/warning messages, give users a hint on how they can see full details.
        """
        return (
            "To view full validation lists, run: "
            "divbase-cli dimensions validate-metadata-file <TSV_FILE_NAME> --project <PROJECT_NAME> "
            "--untruncated."
        )

    def _collect_statistics(self, df: pd.DataFrame) -> ValidationStats:
        """Collect statistics about the TSV file."""

        stats = ValidationStats()
        stats.total_columns = len(df.columns)
        stats.user_defined_columns = len(df.columns) - 1  # Exclude Sample_ID

        project_samples = self.project_samples if self.project_samples is not None else set()
        matching_samples = self.tsv_samples & project_samples
        stats.samples_in_tsv = len(self.tsv_samples)
        stats.samples_matching_project = len(matching_samples)
        stats.total_project_samples = len(project_samples)

        stats.numeric_column_count = len(self.result.numeric_columns)
        stats.string_column_count = len(self.result.string_columns)
        stats.mixed_type_column_count = len(self.result.mixed_type_columns)

        # If has_multi_values is True: at least one cell in the DataFrame contains a Python list (multi-value cell).
        stats.has_multi_values = False
        for col in df.columns:
            for val in df[col].dropna():
                if isinstance(val, list):
                    stats.has_multi_values = True
                    break
            if stats.has_multi_values:
                break

        stats.empty_cells_per_column = {}
        for col in df.columns:
            if col == "Sample_ID":
                continue
            empty_count = 0
            for val in df[col]:
                if val is None or (isinstance(val, str) and val == ""):
                    empty_count += 1
                    continue
                if isinstance(val, list):
                    continue
                if pd.isna(val):
                    empty_count += 1
            if empty_count > 0:
                stats.empty_cells_per_column[col] = empty_count
