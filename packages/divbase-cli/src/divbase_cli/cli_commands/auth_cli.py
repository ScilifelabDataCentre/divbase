"""
CLI subcommand for managing user auth with DivBase server.
"""

from pathlib import Path

import typer
from typing_extensions import Annotated

from divbase_cli.cli_commands.user_config_cli import CONFIG_FILE_OPTION
from divbase_cli.cli_commands.version_cli import PROJECT_NAME_OPTION
from divbase_cli.config_resolver import resolve_project
from divbase_cli.user_auth import login_to_divbase, logout_of_divbase

auth_app = typer.Typer(
    no_args_is_help=True, help="Login/logout of DivBase server. To register, visit https://divbase.scilifelab.se/."
)


@auth_app.command("login")
def login(
    email: str,
    password: Annotated[str, typer.Option(prompt=True, hide_input=True)],
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """
    Log in to the DivBase server.
    """
    project_config = resolve_project(project_name=project, config_path=config_file)
    login_to_divbase(email=email, password=password, divbase_url=project_config.divbase_url)
    print(f"Logged in successfully with email: {email}")


@auth_app.command("logout")
def logout():
    """
    Log out of the DivBase server.
    """
    logout_of_divbase()
    print("Logged out successfully.")
