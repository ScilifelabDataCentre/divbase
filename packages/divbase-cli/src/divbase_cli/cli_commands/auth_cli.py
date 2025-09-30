"""
CLI subcommand for managing user auth with DivBase server.
"""

from pathlib import Path

import typer
from typing_extensions import Annotated

from divbase_cli.cli_commands.user_config_cli import CONFIG_FILE_OPTION
from divbase_cli.config import settings
from divbase_cli.user_auth import login_to_divbase, logout_of_divbase, make_authenticated_request
from divbase_cli.user_config import load_user_config
from divbase_lib.exceptions import AuthenticationError

auth_app = typer.Typer(
    no_args_is_help=True, help="Login/logout of DivBase server. To register, visit https://divbase.scilifelab.se/."
)


@auth_app.command("login")
def login(
    email: str,
    password: Annotated[str, typer.Option(prompt=True, hide_input=True)],
    divbase_url: str = typer.Option(settings.DEFAULT_DIVBASE_API_URL, help="DivBase server URL to connect to."),
    config_file: Path = CONFIG_FILE_OPTION,
):
    """
    Log in to the DivBase server.
    """
    # TODO handle already logged in.
    login_to_divbase(email=email, password=password, divbase_url=divbase_url)
    config = load_user_config(config_file)
    config.set_logged_in_url(divbase_url)

    print(f"Logged in successfully with email: {email}")


@auth_app.command("logout")
def logout():
    """
    Log out of the DivBase server.
    """
    logout_of_divbase()
    print("Logged out successfully.")


@auth_app.command("whoami")
def whoami(
    config_file: Path = CONFIG_FILE_OPTION,
):
    """
    Return information about the currently logged-in user.
    """
    config = load_user_config(config_file)
    logged_in_url = config.logged_in_url
    if not logged_in_url:
        raise AuthenticationError("You are not logged in. Please log in with 'divbase-cli auth login [EMAIL]'.")

    request = make_authenticated_request(
        method="GET",
        divbase_base_url=logged_in_url,
        api_route="v1/auth/whoami",
    )

    if not request.json():
        raise AuthenticationError("You are not logged in. Please log in with 'divbase-cli auth login [EMAIL]'.")

    current_user = request.json()
    print(f"Currently logged in as: {current_user['email']} (Name: {current_user['name']})")
