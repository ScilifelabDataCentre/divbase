"""
To avoid potential problems with circular imports, we can put shared typer args + options (etc...) here

If you have something that is only used in one of the cli subcommands, don't move it here.
"""

import typer

PROJECT_NAME_OPTION = typer.Option(
    None,
    help="Name of the DivBase project, if not provided uses the default in your DivBase config file",
    show_default=False,
)


FORMAT_AS_TSV_OPTION = typer.Option(
    False,
    "--tsv",
    help="If set, will print the output in .TSV format for easier programmatic parsing.",
)
