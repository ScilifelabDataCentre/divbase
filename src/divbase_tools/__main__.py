"""
Manage versions of the entire bucket

This script is interacted with via argparse.

Versioning specified by updating a file in the bucket called ".bucket_versions.yaml".

Versioning done by setting the version according to the current timestamp.
The version file is always kept in the bucket, not stored on disk.

TODOs:
- Check how robust timestamping approach is.
- Add support to delete versions?
- Add support to specify a timestamp for a version?
- Hoist I/O
- How to handle parsing errors/failures upwards to the CLI?
"""

import logging
import sys

import typer

from divbase_tools.cli.user_config_cli import config_app
from divbase_tools.cli.version_cli import version_app

# TODO - check logging config and swap from printing to logging
# and decide how much logging we want to do.
logger = logging.getLogger(__name__)


app = typer.Typer(
    help="""This tool lets you interact with your DivBase project's bucket and do:
        - Upload and download files to/from the bucket.
        - Managed the Version the state of the entire bucket.
    """,
    no_args_is_help=True,
    context_settings={"help_option_names": ["-h", "--help"]},
)


app.add_typer(version_app, name="version")
# app.add_typer(file_app, name="file")
app.add_typer(config_app, name="config")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, handlers=[logging.StreamHandler(sys.stdout)])
    logger.info("Starting divbase_tools CLI application.")
    app()
