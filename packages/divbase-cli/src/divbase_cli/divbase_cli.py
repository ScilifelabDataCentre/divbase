"""
The entry point for the DivBase CLI tool that lets users interact with DivBase project.
"""

import logging
import os
import sys

import typer
from rich import print

from divbase_cli import __version__
from divbase_cli.cli_commands.auth_cli import auth_app
from divbase_cli.cli_commands.dimensions_cli import dimensions_app
from divbase_cli.cli_commands.file_cli import file_app
from divbase_cli.cli_commands.query_cli import query_app
from divbase_cli.cli_commands.task_history_cli import task_history_app
from divbase_cli.cli_commands.user_config_cli import config_app
from divbase_cli.cli_commands.version_cli import version_app
from divbase_cli.cli_config import cli_settings

logger = logging.getLogger(__name__)


app = typer.Typer(
    help="""
    This tool lets you interact with your DivBase project(s) in order to: \n
        - Query the metadata for the VCF files stored in the project. \n
        - Upload/download files to/from the project. \n
        - Version the state of all files in the entire project at a given timestamp.
    """,
    no_args_is_help=True,
    context_settings={"help_option_names": ["-h", "--help"]},
)


def version_callback(show_version: bool) -> None:
    """
    Callback function to display the installed version of the divbase-cli package.
    See https://typer.tiangolo.com/tutorial/options/version/#fix-with-is_eager
    """
    if show_version:
        typer.echo(f"divbase-cli version: {__version__}")
        raise typer.Exit(0)


@app.callback()
def show_installed_version(
    show_version: bool = typer.Option(
        None,
        "--version",
        "-v",
        callback=version_callback,
        is_eager=True,
        help="Show the currently installed version of the divbase-cli package.",
    ),
) -> None:
    pass


app.add_typer(auth_app, name="auth")
app.add_typer(config_app, name="config")
app.add_typer(dimensions_app, name="dimensions")
app.add_typer(file_app, name="files")
app.add_typer(query_app, name="query")
app.add_typer(task_history_app, name="task-history")
app.add_typer(version_app, name="version")


def main():
    # auto-complete mode + logging does not work, as the log message gets included in the auto-complete output.
    # so when running divbase-cli in auto-complete mode we should not log.
    # (when any CLI command is actually run, you are not in_auto_complete_mode)
    in_auto_complete_mode = "_DIVBASE_CLI_COMPLETE" in os.environ
    if cli_settings.LOGGING_ON and not in_auto_complete_mode:
        logging.basicConfig(level=cli_settings.LOG_LEVEL, handlers=[logging.StreamHandler(sys.stderr)])
        logger.info(f"Starting divbase_cli CLI application with logging level: {cli_settings.LOG_LEVEL}")

    # Pretty print any errors for users so they don't see full traceback, unless they have set: DIVBASE_TRACEBACKS_ON=1
    try:
        app()
    except Exception as exc:
        if cli_settings.TRACEBACKS_ON:
            raise
        print(f"[red bold]Error:[/red bold] {str(exc)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
