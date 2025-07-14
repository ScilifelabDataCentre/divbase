"""
A CLI tool to interact with DivBase project(s)/bucket(s).

TODOs:
- Check how robust timestamping approach is.
"""

import logging
import os
import sys
from pathlib import Path

import dotenv
import typer

from divbase_tools.cli_commands.file_cli import file_app
from divbase_tools.cli_commands.query_cli import query_app
from divbase_tools.cli_commands.user_config_cli import config_app
from divbase_tools.cli_commands.version_cli import version_app

logger = logging.getLogger(__name__)


def determine_cli_env() -> None:
    """
    Determine which .env file to use based on enviroment variable passed.
    If not set, will default to taking settings from '.env'.
    """
    ALLOWED_ENVS = ["local", "test", "scilifelab2dev"]
    env_name = os.getenv("DIVBASE_ENV")

    if env_name and env_name not in ALLOWED_ENVS:
        raise ValueError(f"Invalid DIVBASE_ENV value: '{env_name}'. Allowed values are: {', '.join(ALLOWED_ENVS)}")

    if not env_name:
        dotenv_path = Path(".env")
        logger.warning("DIVBASE_ENV not set, using default .env file. ")
    else:
        dotenv_path = Path(f".env.{env_name}")

    if dotenv_path.exists():
        dotenv.load_dotenv(dotenv_path=dotenv_path, override=True)
        logger.info(f"Loaded environment variables from {dotenv_path}")
    else:
        raise FileNotFoundError(
            f"DivBase environment file: '{dotenv_path}' not found. Please create it or set correct ENV variable."
        )


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


def main():
    logging.basicConfig(level=logging.INFO, handlers=[logging.StreamHandler(sys.stdout)])
    logger.info("Starting divbase_tools CLI application.")
    determine_cli_env()
    app()


if __name__ == "__main__":
    main()
