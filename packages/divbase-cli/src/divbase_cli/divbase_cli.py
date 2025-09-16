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

from divbase_cli.cli_commands.dimensions_cli import dimensions_app
from divbase_cli.cli_commands.file_cli import file_app
from divbase_cli.cli_commands.query_cli import query_app
from divbase_cli.cli_commands.user_config_cli import config_app
from divbase_cli.cli_commands.version_cli import version_app

logger = logging.getLogger(__name__)


def load_cli_env_vars() -> None:
    """
    Determine and load in the correct env variables to use for DivBase CLI.
    Env to use controlled by `DIVBASE_ENV` environment variable (can pass in front of any cli command).
    If not set, will default to taking settings from '.env' (for production user).

    A note on the lack of a "test" enviroment. Env variables set in the conftest.py as a fixture.
    The e2e cli tests use the "app" directly (https://typer.tiangolo.com/tutorial/testing/) as recommended.
    They therefore bypass this function as expected.
    """
    ALLOWED_ENVS = ["local", "scilifelab2dev", "scilifelab2prod"]
    env_name = os.getenv("DIVBASE_ENV")
    if not env_name:
        dotenv.load_dotenv(dotenv_path=Path(".env"), override=True)
        logger.warning("DIVBASE_ENV not set, using default .env file. This should probably only be used by an end user")
        return

    if env_name not in ALLOWED_ENVS:
        raise ValueError(f"Invalid DIVBASE_ENV value: '{env_name}'. Allowed values are: {', '.join(ALLOWED_ENVS)}")

    if env_name == "local":
        os.environ["DIVBASE_S3_ACCESS_KEY"] = "minioadmin"
        os.environ["DIVBASE_S3_SECRET_KEY"] = "badpassword"
        logger.info(f"Using hardcoded '{env_name}' CLI environment settings.")
        return

    dotenv_path = Path(f".env.{env_name}")
    if not dotenv_path.exists():
        raise FileNotFoundError(
            f"DivBase environment file: '{dotenv_path}' not found. Please create it or set correct ENV variable."
        )

    dotenv.load_dotenv(dotenv_path=dotenv_path, override=True)
    logger.info(f"Loaded CLI environment variables from {dotenv_path}")


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


def main():
    logging.basicConfig(level=logging.INFO, handlers=[logging.StreamHandler(sys.stdout)])
    logger.info("Starting divbase_cli CLI application.")
    load_cli_env_vars()
    app()


if __name__ == "__main__":
    main()
