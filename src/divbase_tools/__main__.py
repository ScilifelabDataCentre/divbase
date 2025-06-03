"""
A CLI tool to interact with DivBase project(s)/bucket(s).

TODOs:
- Check how robust timestamping approach is.
- Add support to delete versions?
- Add support to specify a timestamp for a version?
"""

import logging
import sys

import typer

from divbase_tools.cli.file_cli import file_app
from divbase_tools.cli.query_cli import query_app
from divbase_tools.cli.user_config_cli import config_app
from divbase_tools.cli.version_cli import version_app

# TODO - check logging config and swap from printing to logging
# and decide how much logging we want to do.
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
app.add_typer(file_app, name="file")
app.add_typer(config_app, name="config")
app.add_typer(query_app, name="query")

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, handlers=[logging.StreamHandler(sys.stdout)])
    logger.debug("Starting divbase_tools CLI application.")
    app()
