"""Collection of utility functions for divbase-cli package that haven't found a better home"""

import csv
import sys

from rich.table import Table


def format_file_size(size_bytes: int | float | None) -> str:
    """
    Converts a file size in bytes to a human-readable format.

    Uses powers of 1000 so KB, MB, GB, TB and not 1024 KiB, MiB, GiB, TiB.
    """
    if size_bytes is None:
        return "N/A"
    if size_bytes == 0:
        return "0 B"
    power = 1000
    n = 0
    power_labels = {0: "", 1: "K", 2: "M", 3: "G", 4: "T"}
    while size_bytes >= power and n < len(power_labels) - 1:
        size_bytes /= power
        n += 1
    return f"{size_bytes:.2f} {power_labels[n]}B"


def print_rich_table_as_tsv(table: Table) -> None:
    """
    Helper function to print a rich Table as a TSV file to standard output.

    This is useful for CLI commands that want to offer both rich table output
    for human users as well as TSV output for programmatic parsing.
    """
    writer = csv.writer(sys.stdout, delimiter="\t")

    headers = [str(col.header) for col in table.columns]
    writer.writerow(headers)

    columns_data = [col._cells for col in table.columns]

    num_rows = len(columns_data[0])
    for row_index in range(num_rows):
        row = [str(columns_data[col_index][row_index]) for col_index in range(len(columns_data))]
        writer.writerow(row)
