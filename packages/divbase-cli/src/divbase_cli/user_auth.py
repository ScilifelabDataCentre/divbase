"""
Manage user authentication with the DivBase server.

This includes login/logout and the getting, storing, using, and refreshing of access + refresh JWTs and Personal Access Tokens (PATs).

User JWTs are stored in the devices's OS keyring. They fallback to a local file with 0600 permissions if not possible.
PATS provided via enviroment variable
"""

import contextlib
import json
import logging
import os
import stat
import time
import warnings
from dataclasses import dataclass
from json import JSONDecodeError
from pathlib import Path

import httpx
import keyring
import stamina
import yaml
from keyring.errors import KeyringError, NoKeyringError, PasswordDeleteError
from pydantic import SecretStr

from divbase_cli import __version__ as cli_version
from divbase_cli.cli_config import cli_settings
from divbase_cli.cli_exceptions import AuthenticationError, DivBaseAPIConnectionError, DivBaseAPIError
from divbase_cli.retries import retry_only_on_retryable_divbase_api_errors
from divbase_cli.user_config import load_user_config
from divbase_lib.api_schemas.auth import LogoutRequest
from divbase_lib.divbase_constants import CLI_VERSION_HEADER_KEY

LOGIN_AGAIN_MESSAGE = "Your session has expired. Please log in again with 'divbase-cli auth login [EMAIL]'."

logger = logging.getLogger(__name__)


@dataclass
class TokenData:
    """Class to hold user JWT data"""

    access_token: SecretStr
    refresh_token: SecretStr
    access_token_expires_at: int
    refresh_token_expires_at: int

    def dump_tokens(self, output_path: Path = cli_settings.TOKENS_PATH) -> None:
        """Dump the user token data to the OS keyring. If fails, fall back to a local file."""
        token_dict = {
            "access_token": self.access_token.get_secret_value(),
            "refresh_token": self.refresh_token.get_secret_value(),
            "access_token_expires_at": self.access_token_expires_at,
            "refresh_token_expires_at": self.refresh_token_expires_at,
        }
        try:
            keyring.set_password(
                service_name=cli_settings.KEYRING_SERVICE,
                username=cli_settings.KEYRING_USERNAME,
                password=json.dumps(token_dict),
            )
            output_path.unlink(missing_ok=True)
            logger.debug("JWTs stored in device keyring successfully.")
            return
        except KeyringError as e:
            logger.debug(f"Keyring JWT storage failed with error: {e}\n falling back to file storage.")

        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as file:
            yaml.safe_dump(token_dict, file, sort_keys=False)
        # only user has read and write permissions for tokens (0600)
        os.chmod(path=output_path, mode=stat.S_IRUSR | stat.S_IWUSR)

    def is_access_token_expired(self) -> bool:
        """Check if the access token is expired"""
        return time.time() >= (self.access_token_expires_at - 5)  # 5 second buffer

    def is_refresh_token_expired(self) -> bool:
        """Check if the refresh token is expired"""
        return time.time() >= (self.refresh_token_expires_at - 300)  # 5 minute buffer


def _delete_stored_jwts(token_path: Path) -> None:
    """Attempt to delete user JWTs from both the keyring and the fallback file."""
    with contextlib.suppress(NoKeyringError, PasswordDeleteError):
        keyring.delete_password(service_name=cli_settings.KEYRING_SERVICE, username=cli_settings.KEYRING_USERNAME)
    token_path.unlink(missing_ok=True)


def check_existing_session(divbase_url: str, config) -> int | None:
    """
    Check if a user is already logged in to DivBase.
    Used to prevent unnecessary multiple logins.

    Returns the refresh token expiry timestamp if logged in (POSIX time), else None.
    """
    if not config.logged_in_url or config.logged_in_url != divbase_url:
        return None

    token_data = load_user_tokens()
    if not token_data or token_data.is_refresh_token_expired():
        return None

    return token_data.refresh_token_expires_at


def _handle_divbase_api_error(response: httpx.Response, http_method: str, url: str) -> None:
    """Handles custom display of a HTTP error response returned by DivBase API."""
    try:
        response_body = response.json()
        error_details = response_body.get("detail", "No error message provided.")
        error_type = response_body.get("type", "unknown")
    except (JSONDecodeError, ValueError):
        # most likely situation for this would be if the reverse proxy returns an error, not the server itself.
        error_details = response.text
        error_type = "unexpected_server_error"

    raise DivBaseAPIError(
        error_details=error_details,
        status_code=response.status_code,
        error_type=error_type,
        http_method=http_method,
        url=url,
    ) from None


@stamina.retry(on=retry_only_on_retryable_divbase_api_errors, attempts=3)
def login_to_divbase(email: str, password: SecretStr, divbase_url: str) -> None:
    """Log in to the DivBase server and return user tokens."""
    login_url = f"{divbase_url}/v1/auth/login"
    try:
        response = httpx.post(
            url=login_url,
            data={
                "grant_type": "password",
                "username": email,  # OAuth2 uses 'username', not 'email'
                "password": password.get_secret_value(),
            },
            headers={
                "Content-Type": "application/x-www-form-urlencoded",
                CLI_VERSION_HEADER_KEY: cli_version,
            },
        )
    except httpx.ConnectError:
        # We don't raise the full error as it contains the password in the stack trace.
        # a user could unknowingly dump this into e.g. a bug report/GitHub issue.
        raise DivBaseAPIConnectionError() from None

    try:
        response.raise_for_status()
    except httpx.HTTPStatusError:
        if response.status_code == 401:
            error_message = response.json().get("detail", "Invalid email or password.")
            raise AuthenticationError(error_message) from None
        _handle_divbase_api_error(response=response, http_method="POST", url=login_url)

    data = response.json()
    token_data = TokenData(
        access_token=SecretStr(data["access_token"]),
        refresh_token=SecretStr(data["refresh_token"]),
        access_token_expires_at=data["access_token_expires_at"],
        refresh_token_expires_at=data["refresh_token_expires_at"],
    )
    token_data.dump_tokens()

    config = load_user_config()
    config.set_login_status(url=divbase_url, email=email)


def logout_of_divbase(token_path: Path = cli_settings.TOKENS_PATH) -> None:
    """
    Log out of the DivBase server.
    We send the refresh token to DivBase to be revoked server-side.
    """
    config = load_user_config()

    # the "if" avoids raising an error on a non logged in user trying to logout
    if config.logged_in_url:
        token_data = load_user_tokens(token_path=token_path)
        if not token_data:
            _delete_stored_jwts(token_path)
            config.set_login_status(url=None, email=None)
            return None

        request_data = LogoutRequest(refresh_token=token_data.refresh_token.get_secret_value())
        # We don't want logout to fail if server is unreachable or gives an error
        # JWTs are stateless so local logout is sufficient.
        try:
            make_unauthenticated_request(
                method="POST",
                divbase_base_url=config.logged_in_url,
                api_route="v1/auth/logout",
                json=request_data.model_dump(),
            )
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

    _delete_stored_jwts(token_path)
    config.set_login_status(url=None, email=None)


def load_user_tokens(token_path: Path = cli_settings.TOKENS_PATH) -> TokenData | None:
    """
    Load user tokens from the OS keyring. Fallback for this is a local file with 0600 permissions.
    Return None if no tokens found or if tokens are malformed/invalid - user logged out.
    """
    try:
        raw = keyring.get_password(service_name=cli_settings.KEYRING_SERVICE, username=cli_settings.KEYRING_USERNAME)
    except KeyringError as e:
        raw = None
        logger.debug(f"Keyring read failed with error: {e}, falling back to file storage.")

    if raw is not None:
        try:
            token_dict = json.loads(raw)
            return TokenData(
                access_token=SecretStr(token_dict["access_token"]),
                refresh_token=SecretStr(token_dict["refresh_token"]),
                access_token_expires_at=token_dict["access_token_expires_at"],
                refresh_token_expires_at=token_dict["refresh_token_expires_at"],
            )
        except (KeyError, json.JSONDecodeError) as e:
            logger.debug(f"Keyring token data malformed: {e}")

    if not token_path.exists():
        return None

    with open(token_path, "r") as file:
        token_dict = yaml.safe_load(file)

    try:
        token_data = TokenData(
            access_token=SecretStr(token_dict["access_token"]),
            refresh_token=SecretStr(token_dict["refresh_token"]),
            access_token_expires_at=token_dict["access_token_expires_at"],
            refresh_token_expires_at=token_dict["refresh_token_expires_at"],
        )
    except KeyError as e:
        logger.debug(f"User tokens file at {token_path} appears to be malformed: {e}")
        return None

    return token_data


@stamina.retry(on=retry_only_on_retryable_divbase_api_errors, attempts=3)
def make_authenticated_request(
    method: str,
    divbase_base_url: str,
    api_route: str,
    token_path: Path = cli_settings.TOKENS_PATH,
    **kwargs,
) -> httpx.Response:
    """Make an authenticated request to the DivBase server, handles refreshing tokens if needed."""
    headers = kwargs.get("headers", {})
    headers[CLI_VERSION_HEADER_KEY] = cli_version

    token_data = load_user_tokens(token_path=token_path)
    if not token_data:
        if not cli_settings.DIVBASE_API_PAT:
            raise AuthenticationError(LOGIN_AGAIN_MESSAGE)
        logger.info("Using personal access token in request to DivBase API.")
        headers["Authorization"] = f"Bearer {cli_settings.DIVBASE_API_PAT.get_secret_value()}"
    else:
        if token_data.is_access_token_expired():
            if token_data.is_refresh_token_expired():
                # Prevents user getting warning about being already logged in when they try to log in again
                config = load_user_config()
                config.set_login_status(url=None, email=None)
                raise AuthenticationError(LOGIN_AGAIN_MESSAGE)
            else:
                token_data = _refresh_access_token(token_data=token_data, divbase_base_url=divbase_base_url)
        headers["Authorization"] = f"Bearer {token_data.access_token.get_secret_value()}"

    kwargs["headers"] = headers
    url = f"{divbase_base_url}/{api_route.lstrip('/')}"

    try:
        response = httpx.request(method=method, url=url, **kwargs)
    except httpx.HTTPError as e:
        raise DivBaseAPIConnectionError() from e

    try:
        response.raise_for_status()
    except httpx.HTTPStatusError:
        _handle_divbase_api_error(response=response, http_method=method, url=url)

    return response


@stamina.retry(on=retry_only_on_retryable_divbase_api_errors, attempts=3)
def make_unauthenticated_request(
    method: str,
    divbase_base_url: str,
    api_route: str,
    **kwargs,
) -> httpx.Response:
    """
    Make an unauthenticated request to the DivBase server.
    Used for those few endpoints that require no authentication, even if the user is logged in.
    E.g. the announcements endpoint.
    """
    url = f"{divbase_base_url}/{api_route.lstrip('/')}"

    headers = kwargs.get("headers", {})
    headers[CLI_VERSION_HEADER_KEY] = cli_version
    kwargs["headers"] = headers

    try:
        response = httpx.request(method=method, url=url, **kwargs)
    except httpx.HTTPError as e:
        raise DivBaseAPIConnectionError() from e

    try:
        response.raise_for_status()
    except httpx.HTTPStatusError:
        _handle_divbase_api_error(response=response, http_method=method, url=url)

    return response


def _refresh_access_token(token_data: TokenData, divbase_base_url: str) -> TokenData:
    """
    Use the refresh token to get a new access token and update the token file.

    Returns the new TokenData object which can be used immediately in a new request.
    NOTE: We do not need retry logic inside this function as the calling function has it.
    """
    refresh_url = f"{divbase_base_url}/v1/auth/refresh"
    try:
        response = httpx.post(
            url=refresh_url,
            json={"refresh_token": token_data.refresh_token.get_secret_value()},
            headers={CLI_VERSION_HEADER_KEY: cli_version},
        )
    except httpx.HTTPError as e:
        raise DivBaseAPIConnectionError() from e

    try:
        response.raise_for_status()
    except httpx.HTTPStatusError:
        # Possible if e.g. token revoked on server side.
        if response.status_code == 401:
            # Clear logged in status in user config as tokens no longer valid.
            # Prevents user getting warning about being already logged in when they try to log in again.
            config = load_user_config()
            config.set_login_status(url=None, email=None)
            raise AuthenticationError(LOGIN_AGAIN_MESSAGE) from None

        _handle_divbase_api_error(response=response, http_method="POST", url=refresh_url)

    data = response.json()

    new_token_data = TokenData(
        access_token=SecretStr(data["access_token"]),
        refresh_token=token_data.refresh_token,  # (refresh_token is still a SecretStr)
        access_token_expires_at=data["expires_at"],
        refresh_token_expires_at=token_data.refresh_token_expires_at,
    )
    new_token_data.dump_tokens()
    return new_token_data
