"""
CLI subcommand for managing user auth with DivBase server.
"""

import logging
import time
from datetime import datetime

import typer
from pydantic import SecretStr
from rich import print

from divbase_cli.cli_config import cli_settings
from divbase_cli.cli_exceptions import DivBaseAPIConnectionError, DivBaseAPIError
from divbase_cli.config_resolver import resolve_url_for_non_project_specific_commands
from divbase_cli.services.announcements import get_and_display_announcements
from divbase_cli.user_auth import (
    PERSONAL_ACCESS_TOKEN_EXPIRED_MESSAGE,
    PATData,
    check_existing_session,
    delete_stored_pat,
    load_stored_user_pat,
    login_to_divbase,
    logout_of_divbase,
    make_authenticated_request,
)
from divbase_cli.user_config import load_user_config
from divbase_lib.divbase_constants import PAT_TOKEN_PREFIX

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

    Removes your locally stored credentials and invalidates the server-side refresh token.
    """
    logout_of_divbase()
    print("Logged out successfully.")


@auth_app.command("add-pat")
def add_pat(
    name: str = typer.Argument(..., help="Name of the personal access token (PAT). E.g. 'work-laptop-pat'"),
    expires_unix_timestamp: int | None = typer.Option(
        None,
        "--expires",
        "-e",
        help="When the personal access token (PAT) expires (if it does), as a unix timestamp.",
        min=1,
    ),
    overwrite_existing: bool = typer.Option(
        False, "--overwrite-existing", "-o", help="Overwrite the existing stored PAT if one already exists."
    ),
):
    """
    Add a personal access token (PAT) to your device for authentication with DivBase.

    PATs are used for authenticating with the DivBase API instead of using your password, and are recommended for use in scripts and pipelines.
    See https://scilifelabdatacentre.github.io/divbase/user-guides/using-divbase-programmatically/ for more details on how to create a personal access token.
    You can only store one personal access token at a time.
    """
    existing = load_stored_user_pat()
    if existing and not overwrite_existing:
        print(
            f"A personal access token named '{existing.name}' is already stored on your device \n"
            f"It is due to expire on: {existing.pat_expiry_formatted()} \n"
            "Append the flag --overwrite-existing (-o) to this command if you want to replace it. \n"
            "You can only store one PAT at a time."
        )
        raise typer.Exit(code=1)

    if expires_unix_timestamp and expires_unix_timestamp < time.time():
        print("The expiry time you entered is in the past. Please enter a valid expiry time for the PAT.")
        raise typer.Exit(code=1)

    pat = SecretStr(
        typer.prompt("please paste your personal access token now", hide_input=True, confirmation_prompt=False)
    )
    if not pat.get_secret_value().startswith(PAT_TOKEN_PREFIX):
        print(
            "It looks like the token you entered is not a valid personal access token. "
            f"Please make sure you copied the entire token, including the '{PAT_TOKEN_PREFIX}' prefix."
        )
        raise typer.Exit(code=1)

    pat_data = PATData(name=name, pat=pat, pat_expires_at=expires_unix_timestamp)
    pat_data.dump_pat_data()
    print(f"Personal access token '{name}' stored successfully.")


@auth_app.command("rm-pat")
def rm_pat():
    """
    Remove the stored personal access token from your device.

    If no such token exists the command will still succeed.
    """
    delete_stored_pat()
    print("If a personal access token was stored, it has now been removed.")


@auth_app.command("pat-info")
def pat_info():
    """
    Display information about the personal access token stored on this device, if any.

    The token will not be displayed for security reasons, but the name and expiry time will be shown.
    """
    pat_data = load_stored_user_pat()
    if not pat_data:
        print("No personal access token stored on this device.")
        return

    print(f"Name: {pat_data.name}")
    print(f"Expires: {pat_data.pat_expiry_formatted()}")

    if pat_data.is_pat_expired():
        print("[red bold]Warning [/red bold]")
        print(PERSONAL_ACCESS_TOKEN_EXPIRED_MESSAGE)


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
