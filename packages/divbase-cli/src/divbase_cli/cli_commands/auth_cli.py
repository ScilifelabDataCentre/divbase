"""
CLI subcommand for managing user auth with DivBase server.
"""

import logging
from datetime import datetime

import typer
from pydantic import SecretStr
from rich import print

from divbase_cli.cli_config import cli_settings
from divbase_cli.cli_exceptions import DivBaseAPIConnectionError, DivBaseAPIError
from divbase_cli.config_resolver import resolve_url_for_non_project_specific_commands
from divbase_cli.services.announcements import get_and_display_announcements
from divbase_cli.user_auth import (
    check_existing_session,
    login_to_divbase,
    logout_of_divbase,
    make_authenticated_request,
)
from divbase_cli.user_config import load_user_config

logger = logging.getLogger(__name__)

divbase_home_url = cli_settings.DIVBASE_API_URL.replace("/api", "")

auth_app = typer.Typer(
    no_args_is_help=True,
    help=f"Login/logout of DivBase server. To register, visit {divbase_home_url}.",
)


@auth_app.command("login")
def login(
    email: str,
    divbase_url: str = typer.Option(cli_settings.DIVBASE_API_URL, help="DivBase server URL to connect to."),
    force: bool = typer.Option(False, "--force", "-f", help="Force login again even if already logged in"),
):
    """
    Log in to the DivBase server.

    You'll be prompted for your password after running the command.
    """
    password: str = typer.prompt("please enter your password", hide_input=True, confirmation_prompt=False)
    secret_password = SecretStr(password)
    del password  # avoid user passwords showing up in error messages etc...

    config = load_user_config()

    if not force:
        session_expires_at = check_existing_session(divbase_url=divbase_url, config=config)
        if session_expires_at:
            print(f"Already logged in to {divbase_url} with email: {config.logged_in_email}.")
            print(f"Session expires: {datetime.fromtimestamp(session_expires_at)}")

            if not typer.confirm("Do you want to login again? This will replace your current session."):
                print("Login cancelled.")
                return

    login_to_divbase(email=email, password=secret_password, divbase_url=divbase_url)

    try:
        get_and_display_announcements(divbase_base_url=divbase_url)
    except (DivBaseAPIError, DivBaseAPIConnectionError):
        # lets not fail the login process if announcements are not working.
        logger.info("Failed to get announcements after login. Error was: ", exc_info=True)

    print(f"Logged in successfully as: {email}")


@auth_app.command("logout")
def logout():
    """
    Log out of the DivBase server.
    """
    logout_of_divbase()
    print("Logged out successfully.")


@auth_app.command("whoami")
def whoami():
    """
    Return information about the currently logged-in user.
    """
    divbase_url = resolve_url_for_non_project_specific_commands()

    request = make_authenticated_request(
        method="GET",
        divbase_base_url=divbase_url,
        api_route="v1/auth/whoami",
    )

    current_user = request.json()
    print(f"Currently logged in as: {current_user['email']} (Name: {current_user['name']})")
