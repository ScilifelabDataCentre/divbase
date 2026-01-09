"""
Manage user authentication with the DivBase server.

This includes login/logout and the getting, storing, using, and refreshing of access + refresh tokens
"""

import time
import warnings
from dataclasses import dataclass
from pathlib import Path

import httpx
import yaml
from pydantic import SecretStr

from divbase_cli.cli_config import cli_settings
from divbase_cli.cli_exceptions import AuthenticationError, DivBaseAPIConnectionError, DivBaseAPIError
from divbase_cli.user_config import load_user_config
from divbase_lib.api_schemas.auth import LogoutRequest

LOGIN_AGAIN_MESSAGE = "Your session has expired. Please log in again with 'divbase-cli auth login [EMAIL]'."


@dataclass
class TokenData:
    """
    Class to hold user token information.
    """

    access_token: SecretStr
    refresh_token: SecretStr
    access_token_expires_at: int
    refresh_token_expires_at: int

    def dump_tokens(self, output_path: Path = cli_settings.TOKENS_PATH) -> None:
        """Dump the user token data to the specified output path"""
        output_path.parent.mkdir(parents=True, exist_ok=True)

        token_dict = {
            "access_token": self.access_token.get_secret_value(),
            "refresh_token": self.refresh_token.get_secret_value(),
            "access_token_expires_at": self.access_token_expires_at,
            "refresh_token_expires_at": self.refresh_token_expires_at,
        }
        with open(output_path, "w") as file:
            yaml.safe_dump(token_dict, file, sort_keys=False)

    def is_access_token_expired(self) -> bool:
        """Check if the access token is expired"""
        return time.time() >= (self.access_token_expires_at - 5)  # 5 second buffer

    def is_refresh_token_expired(self) -> bool:
        """Check if the refresh token is expired"""
        return time.time() >= (self.refresh_token_expires_at - 300)  # 5 minute buffer


def check_existing_session(divbase_url: str, config) -> int | None:
    """
    Check if a user is already logged in to DivBase.
    Used to prevent unnecessary multiple logins.

    Returns the refresh token expiry timestamp if logged in (POSIX time), else None.
    """
    if not config.logged_in_url or config.logged_in_url != divbase_url:
        return None

    try:
        token_data = load_user_tokens()
    except (AuthenticationError, KeyError):
        # e.g. if user manually deleted or modded token file
        return None

    if token_data.is_refresh_token_expired():
        return None

    return token_data.refresh_token_expires_at


def login_to_divbase(email: str, password: SecretStr, divbase_url: str, config_path: Path) -> None:
    """
    Log in to the DivBase server and return user tokens.
    """
    try:
        response = httpx.post(
            f"{divbase_url}/v1/auth/login",
            data={
                "grant_type": "password",
                "username": email,  # OAuth2 uses 'username', not 'email'
                "password": password.get_secret_value(),
            },
            headers={"Content-Type": "application/x-www-form-urlencoded"},
        )
    except httpx.ConnectError:
        # We don't raise the full error as it contains the password in the stack trace.
        # a user could unknowingly dump this into e.g. a bug report/GitHub issue.
        raise DivBaseAPIConnectionError() from None

    if response.status_code == 401:
        error_message = response.json().get("detail", "Invalid email or password.")
        raise AuthenticationError(error_message)

    response.raise_for_status()

    data = response.json()
    token_data = TokenData(
        access_token=SecretStr(data["access_token"]),
        refresh_token=SecretStr(data["refresh_token"]),
        access_token_expires_at=data["access_token_expires_at"],
        refresh_token_expires_at=data["refresh_token_expires_at"],
    )
    token_data.dump_tokens()

    config = load_user_config(config_path)
    config.set_login_status(url=divbase_url, email=email)


def logout_of_divbase(
    token_path: Path = cli_settings.TOKENS_PATH, config_path: Path = cli_settings.CONFIG_PATH
) -> None:
    """
    Log out of the DivBase server.
    We send the refresh token to DivBase to be revoked server-side.
    """
    config = load_user_config(config_path)

    # the "if" avoids raising an error on a non logged in user trying to logout
    if config.logged_in_url:
        token_data = load_user_tokens(token_path=token_path)
        request_data = LogoutRequest(refresh_token=token_data.refresh_token.get_secret_value())

        # We don't want logout to fail if server is unreachable or gives an error
        # JWTs are stateless so local logout is sufficient.
        try:
            make_authenticated_request(
                method="POST",
                divbase_base_url=config.logged_in_url,
                api_route="v1/auth/logout",
                json=request_data.model_dump(),
            )
        except AuthenticationError:
            # Tokens already expired/invalid/revoked etc...
            pass
        except DivBaseAPIConnectionError as e:
            warnings.warn(
                f"Could not connect to DivBase server to log out fully: '{e}'.\n\n"
                "Continuing local logout.\n"
                "You do not need to do anything, but if you see this message often, please let us know.",
                stacklevel=2,
                category=UserWarning,
            )
        except DivBaseAPIError as e:
            warnings.warn(
                f"Received an error message from DivBase server when attempting to logout:"
                f"'{e.error_message=}'. \n\n"
                "Continuing local logout.\n"
                "You do not need to do anything. If you see this message a lot, please let us know.",
                stacklevel=2,
                category=UserWarning,
            )

    token_path.unlink(missing_ok=True)
    config.set_login_status(url=None, email=None)


def load_user_tokens(token_path: Path = cli_settings.TOKENS_PATH) -> TokenData:
    """
    Load user tokens from the specified path.
    """
    if not token_path.exists():
        raise AuthenticationError(
            f"Your access tokens were not found at {token_path}. Please check you are logged in first."
        )

    with open(token_path, "r") as file:
        token_dict = yaml.safe_load(file)

    return TokenData(
        access_token=SecretStr(token_dict["access_token"]),
        refresh_token=SecretStr(token_dict["refresh_token"]),
        access_token_expires_at=token_dict["access_token_expires_at"],
        refresh_token_expires_at=token_dict["refresh_token_expires_at"],
    )


def make_authenticated_request(
    method: str,
    divbase_base_url: str,
    api_route: str,
    token_path: Path = cli_settings.TOKENS_PATH,
    **kwargs,
) -> httpx.Response:
    """
    Make an authenticated request to the DivBase server, handles refreshing tokens if needed.
    """
    token_data = load_user_tokens(token_path=token_path)

    if token_data.is_access_token_expired():
        if token_data.is_refresh_token_expired():
            # Prevents user getting warning about being already logged in when they try to log in again
            config = load_user_config()
            config.set_login_status(url=None, email=None)
            raise AuthenticationError(LOGIN_AGAIN_MESSAGE)
        else:
            token_data = _refresh_access_token(token_data=token_data, divbase_base_url=divbase_base_url)

    headers = kwargs.get("headers", {})
    headers["Authorization"] = f"Bearer {token_data.access_token.get_secret_value()}"
    kwargs["headers"] = headers

    url = f"{divbase_base_url}/{api_route.lstrip('/')}"

    try:
        response = httpx.request(method, url, **kwargs)
    except httpx.HTTPError as e:
        raise DivBaseAPIConnectionError() from e

    try:
        response.raise_for_status()
    except httpx.HTTPStatusError:
        error_details = response.json().get("detail", "No error details provided")
        error_type = response.json().get("type", "unknown")
        raise DivBaseAPIError(
            error_details=error_details,
            status_code=response.status_code,
            error_type=error_type,
            http_method=method,
            url=url,
        ) from None

    return response


def _refresh_access_token(token_data: TokenData, divbase_base_url: str) -> TokenData:
    """
    Use the refresh token to get a new access token and update the token file.

    Returns the new TokenData object which can be used immediately in a new request.
    """
    response = httpx.post(
        f"{divbase_base_url}/v1/auth/refresh",
        json={
            "refresh_token": token_data.refresh_token.get_secret_value(),
        },
    )

    # Possible if e.g. token revoked on server side.
    if response.status_code == 401:
        # Clear logged in status in user config as tokens no longer valid.
        # Prevents user getting warning about being already logged in when they try to log in again.
        config = load_user_config()
        config.set_login_status(url=None, email=None)
        raise AuthenticationError(LOGIN_AGAIN_MESSAGE)

    response.raise_for_status()
    data = response.json()

    new_token_data = TokenData(
        access_token=SecretStr(data["access_token"]),
        refresh_token=token_data.refresh_token,  # (refresh_token is still a SecretStr)
        access_token_expires_at=data["expires_at"],
        refresh_token_expires_at=token_data.refresh_token_expires_at,
    )
    new_token_data.dump_tokens()
    return new_token_data
