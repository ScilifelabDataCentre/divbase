"""
CLI subcommand for managing user auth with DivBase server.
"""

import logging
import sys
from datetime import datetime

import typer
from pydantic import SecretStr
from rich import print

from divbase_cli.cli_config import cli_settings
from divbase_cli.cli_exceptions import AuthenticationError, DivBaseAPIConnectionError, DivBaseAPIError
from divbase_cli.services.announcements import get_and_display_announcements
from divbase_cli.user_auth import (
    check_existing_session,
    login_to_divbase,
    logout_of_divbase,
    make_authenticated_request,
)
from divbase_cli.user_config import load_user_config

logger = logging.getLogger(__name__)

auth_app = typer.Typer(
    no_args_is_help=True, help="Login/logout of DivBase server. To register, visit https://divbase.scilifelab.se/."
)


@auth_app.command("login")
def login(
    email: str,
    divbase_url: str = typer.Option(cli_settings.DIVBASE_API_URL, help="DivBase server URL to connect to."),
    password_stdin: bool = typer.Option(
        False, "--password-stdin", "-p", help="Provide your DivBase password via standard input (STDIN)."
    ),
    force: bool = typer.Option(False, "--force", "-f", help="Force login again even if already logged in"),
):
    """
    Log in to the DivBase server.

    You'll be prompted for your password after running the command.
    Alternatively you can provide your password via standard input (STDIN).
    """
    if password_stdin:
        if sys.stdin.isatty():
            print(
                "[red bold]Error:[/red bold] '--password-stdin' ('-p') flag was set but nothing was piped to STDIN.",
                file=sys.stderr,
            )
            raise typer.Exit(code=1)
        password = sys.stdin.readline().rstrip("\n")
        if not password:
            print(
                "[red bold]Error:[/red bold] '--password-stdin' ('-p') flag was set but an empty password was read from STDIN.",
                file=sys.stderr,
            )
            raise typer.Exit(code=1)
    else:
        password: str = typer.prompt("please enter your password", hide_input=True, confirmation_prompt=False)
    secret_password = SecretStr(password)
    del password  # avoid user passwords showing up in error messages etc...

    config = load_user_config()

    if not force and not password_stdin:
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
    config = load_user_config()
    if not config.logged_in_url:
        raise AuthenticationError("You are not logged in. Please log in with 'divbase-cli auth login [EMAIL]'.")

    request = make_authenticated_request(
        method="GET",
        divbase_base_url=config.logged_in_url,
        api_route="v1/auth/whoami",
    )

    current_user = request.json()
    print(f"Currently logged in as: {current_user['email']} (Name: {current_user['name']})")
