import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from divbase_api.services.vcf_queries import SampleFileMapping
from divbase_api.worker.crud_dimensions import ProjectVCFDimensionsData
from divbase_lib.api_schemas.queries import SampleFileMappingResult, SampleMetadataQueryTaskResult
from divbase_lib.exceptions import (
    SidecarColumnNotFoundError,
    SidecarInvalidFilterError,
    SidecarMetadataFormatError,
    SidecarNoDataLoadedError,
    SidecarSampleIDError,
)
from divbase_lib.metadata_validator import SharedMetadataValidator, ValidationCategory

logger = logging.getLogger(__name__)


###
### Sample metadata query manager and helper functions
###


@dataclass
class NumericFilterContext:
    """
    Context for parsing numeric filter values in SidecarQueryManager.
    """

    key: str
    filter_string_values: str
    is_negated: bool


def run_sidecar_metadata_query(
    file: Path,
    filter_string: str = None,
    project_id: int = None,
    vcf_dimensions_data: ProjectVCFDimensionsData | None = None,
) -> SampleMetadataQueryTaskResult:
    """
    Run a query on a sidecar metadata TSV file and map samples to VCF files.

    Takes vcf_dimensions_data fetched by an API call in the task layer and extracts the unique sample names.

    """
    project_samples = set()
    if vcf_dimensions_data and vcf_dimensions_data.vcf_files:
        for vcf_entry in vcf_dimensions_data.vcf_files:
            sample_names = vcf_entry.samples
            project_samples.update(sample_names)

    sidecar_manager = SidecarQueryManager(file=file, project_samples=project_samples).run_query(
        filter_string=filter_string
    )
    query_message = sidecar_manager.query_message
    warnings = sidecar_manager.warnings
    unique_sample_ids = sidecar_manager.get_unique_values("Sample_ID")

    logger.info(f"Metadata query returned {len(unique_sample_ids)} unique sample IDs")

    if not vcf_dimensions_data or not vcf_dimensions_data.vcf_files:
        error_msg = f"No VCF dimensions data provided for project {project_id}. "
        error_msg += "Please run 'divbase-cli dimensions update' first."
        raise ValueError(error_msg)

    # sets are useful for uniqueness but they are not ordered. The sample-order in bcftools results files becomes non-derterministic if unique_filenames is a set.
    # This was discovered by the results file checksum e2e tests. Use an itermediate set seen_filenames for uniqueness and unique_filenames for as an ordered list.
    sample_and_filename_subset: list[SampleFileMapping] = []
    unique_filenames: list[str] = []
    seen_filenames: set[str] = set()

    for vcf_entry in vcf_dimensions_data.vcf_files:
        filename = vcf_entry.vcf_file_s3_key
        sample_names = vcf_entry.samples

        for sample_id in sample_names:
            if sample_id in unique_sample_ids:
                sample_and_filename_subset.append(SampleFileMapping(sample_id=sample_id, filename=filename))
                if filename not in seen_filenames:
                    seen_filenames.add(filename)
                    unique_filenames.append(filename)

    return SampleMetadataQueryTaskResult(
        sample_and_filename_subset=[
            SampleFileMappingResult(sample_id=entry.sample_id, filename=entry.filename)
            for entry in sample_and_filename_subset
        ],
        unique_sample_ids=list(unique_sample_ids),
        unique_filenames=unique_filenames,
        query_message=query_message,
        warnings=warnings,
    )


class SidecarQueryManager:
    """
    A class that manages the execution queries on sidecar metadata files.

    Expects a TSV file with a header row and tab-separated values.
    The class provides methods to load the TSV file into a pandas DataFrame, run queries against the data,
    and retrieve unique values from specific columns.

    TODO - consider seperation of concerns: this class currently handles both loading the TSV file and running queries on it.
    TODO - some of the __init__ params are perhaps better as properties?
    """

    def __init__(self, file: Path, project_samples: set[str] | None = None):
        self.file = file
        self.project_samples = project_samples
        self.filter_string = None
        self.df = None
        self.metadata_validator = None
        self.query_result = None
        self.query_message: str = ""
        self.numeric_columns: list[str] = []
        self.string_columns: list[str] = []
        self.mixed_type_columns: list[str] = []
        self.warnings: list[str] = []
        self.load_file()

    def load_file(self) -> "SidecarQueryManager":
        """
        Method that loads the TSV file into a pandas DataFrame. Assumes that the first row is a header row, and that the file is tab-separated.
        Also removes any leading '#' characters from the column names.

        Uses the warning and error category Enums from SharedMetadataValidator logic to raise errors or send warnings to the user.
        """
        try:
            logger.info(f"Loading sidecar metadata file: {self.file}")

            # Note! The SharedMetadataValidator is for checks on the contents of the TSV file. The logic is shared between this class and the client-side ClientSideMetadataTSVValidator.
            # There are several helper methods for the filtering logic in this class, but they are for the query filters and are not related to the validation of the TSV file contents.
            self.metadata_validator = SharedMetadataValidator(
                file_path=self.file,
                project_samples=self.project_samples,
                skip_dimensions_check=(self.project_samples is None),
            )
            result = self.metadata_validator.load_and_validate()

            if result.errors:
                # Note! The order of these errors matters. The first error in the list is the one that is raised, so more critical errors should be placed higher in the order than less critical errors.
                first_encountered_error = result.errors[0]
                if first_encountered_error.category == ValidationCategory.FILE_READ:
                    raise SidecarNoDataLoadedError(file_path=self.file, submethod="load_file")
                elif first_encountered_error.category == ValidationCategory.SAMPLE_ID_COLUMN:
                    raise SidecarColumnNotFoundError(first_encountered_error.message)
                elif first_encountered_error.category == ValidationCategory.SAMPLE_ID_VALUE:
                    raise SidecarSampleIDError(first_encountered_error.message)
                else:
                    raise SidecarMetadataFormatError(first_encountered_error.message)

            if result.warnings:
                self.warnings.extend(
                    w.message
                    for w in result.warnings
                    if w.category in (ValidationCategory.DIMENSIONS, ValidationCategory.FORMAT)
                )

            self.df = result.df
            self.numeric_columns = result.numeric_columns
            self.string_columns = result.string_columns
            self.mixed_type_columns = result.mixed_type_columns

        except (
            SidecarSampleIDError,
            SidecarColumnNotFoundError,
            SidecarInvalidFilterError,
            SidecarMetadataFormatError,
            SidecarNoDataLoadedError,
        ):
            # Let validation errors propagate directly to user with specific error messages
            raise
        except Exception as e:
            # Only wrap unexpected errors (file I/O, pandas errors, etc.)
            raise SidecarNoDataLoadedError(file_path=self.file, submethod="load_file") from e
        return self

    def get_unique_values(self, column: str) -> list:
        """
        Method to fetch unique values from a specific column in the query result. Intended to be invoked on a SidecarQueryManager
        instance after a query has been run with run_query().
        """
        if self.query_result is None:
            raise SidecarColumnNotFoundError("No query result available. Run run_query() first.")

        if column in self.query_result.columns:
            return self.query_result[column].unique().tolist()
        else:
            raise SidecarColumnNotFoundError(f"Column '{column}' not found in query result")

    def run_query(self, filter_string: str = None) -> "SidecarQueryManager":
        """
        Method to run a query against the TSV data loaded by self.load_file(). The filter_string should be a semicolon-separated list of key:value pairs,
        where key is a column name and value is a comma-separated list of filter values.
        For example: "key1:value1,value2;key2:value3,value4".

        The TSV that is loaded into the pandas DataFrame can have both string and numeric columns.
        - String columns are matched to filter string values with OR logic: if ANY value in a cell matches ANY filter value, the row matches.
        - Numeric columns support:
            - Inequalities: ">25", "<=40" (checks if any cell value satisfies the condition)
            - Ranges: "20-40" (checks if any cell value is within the range)
            - Discrete values: "25,30,50" (checks if any cell value matches any filter value)
            - All are combined with OR logic

        Filtering using the ! (NOT) operator:
        - "!" must prefix the filter value, e.g. "key:!value" means that rows with "value" in the "key" column should be excluded.
        - Numeric examples: "Population:!2" (exclude 2), "Age:<30,!25" (less than 30 but not 25), "Weight:!20-40" (exclude range 20-40)
        - String examples: "Area:!North" (exclude North), "Region:East,West,!South" (East or West but not South)
        - NOT conditions are applied with AND logic after positive conditions have been applied: rows must NOT match any negated value

        Filter string values in the query vs. cell values in the TSV:
        - Filter strings are handled per semicolon-separated key-value pair: in "key1:value1,value2;key2:value3,value4"
          "key1:value1,value2" is handled separately from "key2:value3,value4".
        - Filter string values can be comma-separated, e.g. "value1,value2" in "key1:value1,value2" and each filter string value is handled separately.
        - Cells can have multi-values as long as Python list syntax is used in the TSV cell, e.g. [25, 30, 35].
        - Matching of filter string to cell values uses OR logic: if ANY value in a cell matches ANY filter value, the row matches.

        Note:
        Even though multi-value cells now use Python list syntax (e.g., [1, 2, 3]), the validator and query logic still check for semicolons and commas in plain string cells.
        This is  because: semicolons and commas have special meaning in the query filter syntax (semicolon separates key-value pairs, comma separates filter values).
        If a user enters a semicolon or comma in a plain string cell (not a list), it may cause confusion or unexpected query results, as the query parser may split on these characters and thus will not be able to match against them.
        Warnings are issued to help users avoid ambiguous or unintended filter behavior.

        Summary of how different input filter values are handled:
        - If the filter_string is empty, all records are returned.
        - If the filter_string is None, an error is raised.
        - If the filter_string is not empty, the method filters the dataframe based on the provided filter_string.
        - If any of the keys in the filter_string are not found in the dataframe columns, a warning is logged and those conditions are skipped.
        - If none of the filter string values match any cell values in the dataframe, a warning is logged and all records are returned.
        - If the filter_string is invalid, a SidecarInvalidFilterError is raised.

        The method returns the SidecarQueryManager instance with the query_result and query_message. The former is the filtered DataFrame results,
        and the latter is filter_string used for the query.
        """

        if self.df is None:
            raise SidecarNoDataLoadedError(
                file_path=self.file, submethod="run_query", error_details="No data loaded. Call load_file() first."
            )

        if filter_string is not None:
            self.filter_string = filter_string

        if self.filter_string == "":
            logger.warning("Empty filter provided - returning ALL records. This may be a large result set.")
            self.query_result = self.df
            self.query_message = "ALL RECORDS (no filter)"
            return self

        if self.filter_string is None:
            raise SidecarInvalidFilterError("Filter cannot be None. Use an empty string ('') if you want all records.")

        key_values = self.filter_string.split(";")
        filter_conditions = []

        # 1. Parse the input filter string and build a list of boolean conditions to apply to the dataframe
        for key_value in key_values:
            if not key_value.strip():
                continue
            try:
                key, filter_string_values = key_value.split(":", 1)
                key = key.strip()
                filter_string_values = filter_string_values.strip()

                if key not in self.df.columns:
                    warning_msg = f"Column '{key}' not found in the TSV file. Skipping this filter condition."
                    logger.warning(warning_msg)
                    self.warnings.append(warning_msg)
                    continue

                is_numeric = key in self.numeric_columns

                # Check if type consistency and return warnings to users if applicable.
                # If the column is treated as string, check for potential user mistakes (e.g. using numeric filter syntax on a string column that contains numeric-looking values):
                # 1. Warn if the column has mixed types (some values look numeric) and that the column will be treat as string type.
                # 2. Warn if the filter uses numeric syntax on this string column. Do not raise error.
                if not is_numeric:
                    is_mixed = key in self.mixed_type_columns

                    problematic_filter_values = self._detect_numeric_filter_syntax_on_string_column(
                        key, filter_string_values
                    )

                    if is_mixed or problematic_filter_values:
                        warning_lines = [f"Column '{key}':"]
                        if is_mixed:
                            warning_lines.append("      - Contains mixed types (both numeric and non-numeric values).")
                            warning_lines.append("        This column is treated as a string column.")
                        if problematic_filter_values:
                            warning_lines.append(
                                f"      - Your filter contains comparison operators {problematic_filter_values}, "
                                "which are not supported on string columns."
                            )
                            warning_lines.append(
                                "        DivBase comparison operators (>, <, >=, <=) only work on numeric columns."
                            )
                            warning_lines.append(
                                f"        Use exact string matching instead (e.g., '{key}:value1,value2')."
                            )
                        warning_msg = "\n".join(warning_lines)
                        logger.warning(warning_msg)
                        self.warnings.append(warning_msg)

                # Handle numeric filtering: inequalities, ranges, and discrete values (all with OR logic)
                # e.g., "Weight:>25,<30,50" or "Weight:20-40,50,>100"
                # Supports filtering on semicolon-separated values in cells in the TSV: e.g. "25;30;35"
                # Also handles columns that pandas infers as strings but contain numeric values with semicolons (e.g., "1;2;3")
                # Also supports NOT operator with ! prefix: e.g., "Weight:!25" or "Weight:<4,!2"
                if is_numeric:
                    filter_string_values_list = self._split_filter_values(filter_string_values)

                    # Negated values are those that start with "!" in the filter string
                    positive_values, negated_values = self._separate_positive_and_negated_values(
                        filter_values=filter_string_values_list
                    )

                    positive_filter_context = NumericFilterContext(
                        key=key,
                        filter_string_values=filter_string_values,
                        is_negated=False,
                    )

                    inequality_conditions, range_conditions, discrete_values = self._parse_numeric_filter_values(
                        values_to_process=positive_values,
                        context=positive_filter_context,
                    )

                    negated_filter_context = NumericFilterContext(
                        key=key,
                        filter_string_values=filter_string_values,
                        is_negated=True,
                    )
                    negated_inequality_conditions, negated_range_conditions, negated_discrete_values = (
                        self._parse_numeric_filter_values(
                            values_to_process=negated_values,
                            context=negated_filter_context,
                        )
                    )

                    # Combine multiple conditions (inequality, range, discrete values) with OR logic
                    conditions = self._build_condition_list(
                        inequality_conditions=inequality_conditions,
                        range_conditions=range_conditions,
                        discrete_values=discrete_values,
                        key=key,
                    )

                    negated_conditions = self._build_condition_list(
                        inequality_conditions=negated_inequality_conditions,
                        range_conditions=negated_range_conditions,
                        discrete_values=negated_discrete_values,
                        key=key,
                    )

                    if conditions or negated_conditions:
                        # First combine all positive conditions with OR logic. Can be None if there are no positive conditions, e.g. if the filter string only contains negated conditions like "Weight:!20-40"
                        base_condition = self._combine_conditions_with_or(conditions=conditions) if conditions else None
                        # Then apply negated conditions with AND logic: rows must NOT match any negated condition.
                        combined = self._apply_not_conditions(
                            base_condition=base_condition, negated_conditions=negated_conditions
                        )

                        if not combined.any():
                            warning_msg = f"No values in column '{key}' match the filter: {filter_string_values}"
                            logger.warning(warning_msg)
                            self.warnings.append(warning_msg)
                        filter_conditions.append(combined)

                    else:
                        warning_msg = f"No valid numeric values, ranges, or inequalities provided for column '{key}'. Filter condition will not match any rows."
                        logger.warning(warning_msg)
                        self.warnings.append(warning_msg)
                else:
                    # Non-numeric column: handle as discrete string values
                    # Supports NOT operator with ! prefix: e.g., "Area:!North" or "Area:North,!South"
                    # Supports quoted values with commas: e.g., 'Area:"North,South"' matches the literal string
                    filter_string_values_list = self._split_filter_values(filter_string_values)

                    positive_values, negated_values = self._separate_positive_and_negated_values(
                        filter_values=filter_string_values_list
                    )

                    # Build condition
                    if positive_values or negated_values:
                        base_condition = (
                            self._create_string_condition(key=key, target_values=positive_values)
                            if positive_values
                            else None
                        )
                        negated_conditions = (
                            [self._create_string_condition(key=key, target_values=negated_values)]
                            if negated_values
                            else []
                        )

                        condition = self._apply_not_conditions(
                            base_condition=base_condition, negated_conditions=negated_conditions
                        )

                        if not condition.any():
                            warning_msg = (
                                f"No results for the filter {filter_string_values_list} were found in column '{key}'."
                            )
                            logger.warning(warning_msg)
                            self.warnings.append(warning_msg)
                        filter_conditions.append(condition)
            except SidecarInvalidFilterError:
                # Allow specific validation errors (like "contains commas") to propagate unchanged.
                # This preserves detailed error messages for user-facing exceptions.
                raise
            except Exception as e:
                raise SidecarInvalidFilterError(
                    f"Invalid filter format: '{key_value}'. Expected format 'key:value1,value2' or 'key:min-max' for numeric ranges"
                ) from e

        # 2. Apply the final boolean filters on the dataframe
        if filter_conditions:
            combined_condition = pd.Series(True, index=self.df.index)
            # Iteratively combine each condition in filter_conditions to create a final combined condition where each row must satisfy all filter conditions to be included.
            for condition in filter_conditions:
                # The ampersand (&) is pandas syntax for element-wise AND between boolean Series.
                combined_condition = combined_condition & condition

            self.query_result = self.df[combined_condition].copy()
            self.query_message = self.filter_string
        else:
            raise SidecarInvalidFilterError(
                f"Invalid filter conditions: no valid filter conditions could be parsed from '{self.filter_string}'. "
                "Please check your filter keys, value spelling, and syntax. "
                "Expected format: 'Key:Value1,Value2' or 'Key:min-max' for numeric ranges."
            )

        return self

    def _detect_numeric_filter_syntax_on_string_column(self, key: str, filter_string_values: str) -> list[str]:
        """
        Helper method for the filtering logic to detect when a user's filter string contains inequality operators
        (e.g. ">25", ">=10", "<North", "<=abc") on a column that is treated as string.

        Detects any use of inequality operators (>, <, >=, <=) as a prefix, regardless of whether the
        value after the operator is numeric or not. E.g. ">5" and ">North".


        Doesn't flag range-like filter values like "1-2" since these are common string values
        (e.g., hyphenated IDs or names) and will correctly match via string matching.

        Returns a list of the problematic filter values for use in a warning messages.
        """
        problematic_filter_values = []
        values = self._split_filter_values(filter_string_values)
        for filter_value in values:
            filter_value = filter_value.strip().lstrip("!")  # strip negation prefix for checking
            if not filter_value:
                continue
            # Check for inequality operators:
            if re.match(r"^(>=|<=|>|<).+$", filter_value):
                problematic_filter_values.append(filter_value)
        return problematic_filter_values

    def _split_filter_values(self, filter_values_str: str) -> list[str]:
        """
        Split comma-separated filter-value strings.

        Designed to handle cases with filter values that contain commas or other special characters by allowing users to wrap such values in double quotes.
        For example, if a TSV cell contains the literal string: `North, South`, the CLI filter must be wrapped in double quotes so the comma is not
        treated as a value separator:

            divbase-cli query tsv 'Area:"North, South"'

        Additionally, the filter string itself must be wrapped in single quotes to prevent the shell from interpreting the inner double quotes.
        Only double-quote quoting is supported inside filter values; single quotes inside the filter string are treated as literal characters
        (i.e. ``"Area:'North, South'"`` is not supported.).

        The ``!`` NOT operator is preserved as part of the string and is handled later by ``_separate_positive_and_negated_values``.

        Examples:
            "North,South"        -> ["North", "South"]
            '"North,South",East' -> ["North,South", "East"]
            '!"North,South"'    -> ["!North,South"]
        """
        values = []
        current = []
        in_quotes = False

        # Iterate through the filter values string character by character to handle double-quoted strings correctly.
        for char in filter_values_str:
            if char == '"':
                in_quotes = not in_quotes
            elif char == "," and not in_quotes:
                values.append("".join(current).strip())
                current = []
            else:
                current.append(char)

        values.append("".join(current).strip())
        return [v for v in values if v]

    def _get_cell_values(self, cell_value: Any) -> list:
        """Return cell value as a list of values for filtering.

        Checks for list type before pd.isna() because pd.isna() raises
        ValueError on list/array inputs.
        """
        if isinstance(cell_value, list):
            return cell_value
        if pd.isna(cell_value):
            return []
        return [cell_value]

    def _parse_numeric_value(self, value_str: str) -> float | int:
        """Helper method for the filtering logic to parse a string value to int or float. To be used when other checks have already confirmed that the value can be parsed as numeric."""
        return float(value_str) if "." in value_str else int(value_str)

    def _create_inequality_condition(self, key: str, operator: str, threshold: float) -> pd.Series:
        """
        Helper method for the filtering logic to create a condition for inequality filtering on a column.
        Uses a named nested function instead of a lambda to improve readability to defined the
        logic that will be applied to the Pandas dataframe.
        """

        def check_inequality(cell_value):
            cell_values = self._get_cell_values(cell_value)
            if not cell_values:
                return False
            for val in cell_values:
                try:
                    val_num = float(val)
                    if (
                        (operator == ">" and val_num > threshold)
                        or (operator == ">=" and val_num >= threshold)
                        or (operator == "<" and val_num < threshold)
                        or (operator == "<=" and val_num <= threshold)
                    ):
                        return True
                except ValueError:
                    continue
            return False

        return self.df[key].apply(check_inequality)

    def _create_range_condition(self, key: str, min_val: float, max_val: float) -> pd.Series:
        """
        Helper method for the filtering logic to create a condition for range filtering on a column.
        Uses a named nested function instead of a lambda to improve readability to define the
        logic that will be applied to the Pandas dataframe.
        """

        def check_range(cell_value):
            cell_values = self._get_cell_values(cell_value)
            if not cell_values:
                return False
            for val in cell_values:
                try:
                    val_num = float(val)
                    if min_val <= val_num <= max_val:
                        return True
                except ValueError:
                    continue
            return False

        return self.df[key].apply(check_range)

    def _create_discrete_numeric_condition(self, key: str, target_values: list[float | int]) -> pd.Series:
        """
        Helper method for the filtering logic to create a condition for discrete numeric value filtering on a column.
        Uses a named nested function instead of a lambda to improve readability to define the
        logic that will be applied to the Pandas dataframe.
        """

        def check_discrete(cell_value):
            cell_values = self._get_cell_values(cell_value)
            if not cell_values:
                return False
            for val in cell_values:
                try:
                    val_num = float(val)
                    if val_num in target_values:
                        return True
                except ValueError:
                    continue
            return False

        return self.df[key].apply(check_discrete)

    def _create_string_condition(self, key: str, target_values: list[str]) -> pd.Series:
        """
        Helper method for the filtering logic to create a condition for string value filtering on a column.
        Uses a named nested function instead of a lambda to improve readability to define the
        logic that will be applied to the Pandas dataframe.
        """

        def check_string(cell_value):
            cell_values = self._get_cell_values(cell_value)
            return any(str(val) in target_values for val in cell_values)

        return self.df[key].apply(check_string)

    def _combine_conditions_with_or(self, conditions: list[pd.Series]) -> pd.Series:
        """
        Helper method for the filtering logic to combine multiple Pandas boolean Series with OR logic.
        Returns a single boolean Series that is True if any of the input conditions is True for each row.

        The resulting Series is used at the end of self.run_query() to filter the DataFrame values.
        """
        if not conditions:
            return pd.Series(False, index=self.df.index)
        combined = conditions[0]
        for cond in conditions[1:]:
            # The bar (|) is pandas syntax for element-wise OR between boolean Series.
            combined = combined | cond
        return combined

    def _build_condition_list(
        self,
        inequality_conditions: list[pd.Series],
        range_conditions: list[pd.Series],
        discrete_values: list[float | int],
        key: str,
    ) -> list[pd.Series]:
        """
        Helper method for the filtering logic to build a list of conditions from inequality, range, and discrete value filters.
        """
        conditions = []

        if inequality_conditions:
            conditions.append(self._combine_conditions_with_or(conditions=inequality_conditions))

        if range_conditions:
            conditions.append(self._combine_conditions_with_or(conditions=range_conditions))

        if discrete_values:
            discrete_condition = self._create_discrete_numeric_condition(key=key, target_values=discrete_values)
            conditions.append(discrete_condition)

        return conditions

    def _parse_numeric_filter_values(
        self, values_to_process: list[str], context: NumericFilterContext
    ) -> tuple[list[pd.Series], list[pd.Series], list[float | int]]:
        """
        Helper method for the filtering logic to identify if a numeric filter values is an inequality, range, or discrete value and process it accordingly.

        The context dataclass is intended to keep the kwargs manageable when passing positive and negative values back-to-back:
            - key: Column name being filtered
            - filter_string_values: Original filter string (for error messages)
            - is_negated: Whether these are negated (NOT) conditions
        """
        key = context.key
        filter_string_values = context.filter_string_values
        is_negated = context.is_negated

        inequality_conditions = []
        range_conditions = []
        discrete_values = []

        for filter_string_value in values_to_process:
            # Check for common mistakes: =< or => instead of <= or >=
            if re.match(r"^=<-?\d+\.?\d*$", filter_string_value) or re.match(r"^=>-?\d+\.?\d*$", filter_string_value):
                raise SidecarInvalidFilterError(
                    f"Invalid operator format '{filter_string_value[:2]}' in filter '{key}:{filter_string_values}'."
                    f"Use standard operators: '<=' (not '=<') or '>=' (not '=>')"
                )

            # Check if it's an inequality (e.g., ">25", "<=40", "<-5")
            inequality_match = re.match(r"^(>=|<=|>|<)(-?\d+\.?\d*)$", filter_string_value)
            if inequality_match:
                operator = inequality_match.group(1)
                threshold = float(inequality_match.group(2))
                condition = self._create_inequality_condition(key, operator, threshold)
                inequality_conditions.append(condition)
                prefix = "NOT " if is_negated else ""
                logger.debug(
                    f"Applied {'negated ' if is_negated else ''}inequality filter on '{key}': {prefix}{operator} {threshold}"
                )
                continue

            # Check if it's a range (e.g., "20-40", "-100--50", "10-20")
            range_match = re.match(r"^(-?\d+\.?\d*)-(-?\d+\.?\d*)$", filter_string_value)
            if range_match:
                min_val = float(range_match.group(1))
                max_val = float(range_match.group(2))
                condition = self._create_range_condition(key, min_val, max_val)
                range_conditions.append(condition)
                prefix = "NOT " if is_negated else ""
                logger.debug(
                    f"Applied {'negated ' if is_negated else ''}range filter on '{key}': {prefix}{min_val} to {max_val}"
                )
                continue

            # Otherwise, treat as discrete value
            try:
                numeric_value = float(filter_string_value) if "." in filter_string_value else int(filter_string_value)
                discrete_values.append(numeric_value)
            except ValueError:
                logger.warning(
                    f"Cannot convert '{filter_string_value}' to numeric for column '{key}'. Skipping this value."
                )

        return inequality_conditions, range_conditions, discrete_values

    def _separate_positive_and_negated_values(self, filter_values: list[str]) -> tuple[list[str], list[str]]:
        """
        Helper method for the filtering logic to separate filter values into positive and negated lists.
        Values prefixed with '!' are negated (NOT conditions).
        """
        positive_values = []
        negated_values = []

        for value in filter_values:
            value = value.strip()
            if value.startswith("!"):
                negated_values.append(value[1:].strip())
            else:
                positive_values.append(value)

        return positive_values, negated_values

    def _apply_not_conditions(self, base_condition: pd.Series | None, negated_conditions: list[pd.Series]) -> pd.Series:
        """
        Helper method for the filtering logic to apply NOT conditions to a base condition. The base condition contains positive filters combined with OR, or None if there were only negations
        in the input filter string from the CLI. Returns a combined condition where rows must match base_condition AND NOT match any negated condition
        """
        if base_condition is None:
            # If only negated conditions (no positive conditions), start with all True
            combined = pd.Series(True, index=self.df.index)
        else:
            combined = base_condition

        # Apply negated conditions (must NOT match any negated condition)
        for negated_condition in negated_conditions:
            combined = combined & ~negated_condition

        return combined
