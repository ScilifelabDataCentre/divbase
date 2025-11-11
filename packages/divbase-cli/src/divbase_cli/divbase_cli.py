"""
The entry point for the DivBase CLI tool that lets users interact with DivBase project.
"""

import logging
import sys

import typer

from divbase_cli.cli_commands.auth_cli import auth_app
from divbase_cli.cli_commands.dimensions_cli import dimensions_app
from divbase_cli.cli_commands.file_cli import file_app
from divbase_cli.cli_commands.query_cli import query_app
from divbase_cli.cli_commands.user_config_cli import config_app
from divbase_cli.cli_commands.version_cli import version_app

logger = logging.getLogger(__name__)


app = typer.Typer(
    help="""
    This tool lets you interact with your DivBase project(s) bucket(s) in order to: \n
        - Query the metadata for the VCF files stored in the bucket. \n
        - Upload/download files to/from the bucket. \n
        - Version the state of all files in the entire bucket at a given timestamp.
    """,
    no_args_is_help=True,
    context_settings={"help_option_names": ["-h", "--help"]},
)


app.add_typer(version_app, name="version")
app.add_typer(file_app, name="files")
app.add_typer(config_app, name="config")
app.add_typer(query_app, name="query")
app.add_typer(dimensions_app, name="dimensions")
app.add_typer(auth_app, name="auth")


def main():
    logging.basicConfig(level=logging.INFO, handlers=[logging.StreamHandler(sys.stdout)])
    logger.info("Starting divbase_cli CLI application.")
    app()


if __name__ == "__main__":
    main()
