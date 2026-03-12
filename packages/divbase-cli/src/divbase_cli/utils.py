"""Collection of utility functions for divbase-cli package that haven't found a better home"""

import csv
import sys

from rich.table import Table


def print_rich_table_as_tsv(table: Table) -> None:
    """
    Helper function to print a rich Table as a TSV file to standard output.

    This is useful for CLI commands that want to offer both rich table output
    for human users as well as TSV output for programmatic parsing.

    NOTE: This function expects all table rows to be of same length (you can have None values in cells).
    """
    writer = csv.writer(sys.stdout, delimiter="\t")

    headers = [str(col.header) for col in table.columns]
    writer.writerow(headers)

    columns_data = [col._cells for col in table.columns]

    num_rows = min(len(col) for col in columns_data)
    for row_index in range(num_rows):
        row = [str(columns_data[col_index][row_index]) for col_index in range(len(columns_data))]
        writer.writerow(row)
