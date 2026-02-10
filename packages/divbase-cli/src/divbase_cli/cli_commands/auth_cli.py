"""
CLI subcommand for managing user auth with DivBase server.
"""

import logging
from datetime import datetime

import typer
from pydantic import SecretStr
from typing_extensions import Annotated

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
    password: Annotated[str, typer.Option(prompt=True, hide_input=True)],
    divbase_url: str = typer.Option(cli_settings.DIVBASE_API_URL, help="DivBase server URL to connect to."),
    force: bool = typer.Option(False, "--force", "-f", help="Force login again even if already logged in"),
):
    """
    Log in to the DivBase server.

    TODO - think abit more about already logged in validation and UX.
    One thing to consider would be use case of very close to refresh token expiry, that could be bad UX.
    (But that is dependent on whether we will allow renewal of refresh tokens...)
    """
    secret_password = SecretStr(password)
    del password  # avoid user passwords showing up in error messages etc...

    config = load_user_config()

    if not force:
        session_expires_at = check_existing_session(divbase_url=divbase_url, config=config)
        if session_expires_at:
            print(f"Already logged in to {divbase_url}")
            print(f"Session expires: {datetime.fromtimestamp(session_expires_at)}")

            if not typer.confirm("Do you want to login again? This will replace your current session."):
                print("Login cancelled.")
                return

    login_to_divbase(email=email, password=secret_password, divbase_url=divbase_url)

    try:
        get_and_display_announcements(divbase_base_url=divbase_url)
    except (DivBaseAPIError, DivBaseAPIConnectionError):
        # lets not fail the login process if annoucements are not working.
        logger.error("Failed to get announcements after login. Error was: ", exc_info=True)

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
    logged_in_url = config.logged_in_url

    # TODO - move logged in check to the make_authenticated_request function?
    if not logged_in_url:
        raise AuthenticationError("You are not logged in. Please log in with 'divbase-cli auth login [EMAIL]'.")

    request = make_authenticated_request(
        method="GET",
        divbase_base_url=logged_in_url,
        api_route="v1/auth/whoami",
    )

    current_user = request.json()
    print(f"Currently logged in as: {current_user['email']} (Name: {current_user['name']})")
