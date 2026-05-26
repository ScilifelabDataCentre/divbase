"""
To avoid potential problems with circular imports, we can put shared typer args + options (etc...) here

If you have something that is only used in one of the cli subcommands, don't move it here.
"""

import typer

DOWNLOAD_DIR_OPTION = typer.Option(
    None,
    "--download-dir",
    "-d",
    help="""Directory to download the files to.
        If not provided, defaults to what you specified in your user config.
        If also not specified in your user config, downloads to the current directory.
        You can also specify "." to download to the current directory.""",
)

PROJECT_NAME_OPTION = typer.Option(
    None,
    "--project",
    "-p",
    help="Name of the DivBase project, if not provided uses the default in your DivBase config file",
    show_default=False,
)


FORMAT_AS_TSV_OPTION = typer.Option(
    False,
    "--tsv",
    "-t",
    help="If set, will print the output in .TSV format for easier programmatic parsing.",
)
