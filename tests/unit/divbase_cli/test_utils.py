"""
Unit tests for divbase_cli.utils
"""

from rich.table import Table

from divbase_cli.utils import print_rich_table_as_tsv


def test_print_rich_table_as_tsv_keeps_each_row_on_one_line(capsys):
    """
    Cells with newlines or tabs (e.g. multi-line error messages in task history) should not split a row
    over several lines/columns, so the output can be parsed with line based tools like cut and awk.
    """
    table = Table()
    table.add_column("ID")
    table.add_column("Result")
    table.add_row("1", "First line of error.\nSecond line\twith a tab.\r\nThird line.")
    table.add_row("2", "Success")

    print_rich_table_as_tsv(table=table)

    lines = capsys.readouterr().out.splitlines()
    assert lines == [
        "ID\tResult",
        "1\tFirst line of error. Second line with a tab. Third line.",
        "2\tSuccess",
    ]
