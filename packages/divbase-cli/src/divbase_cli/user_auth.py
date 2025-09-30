"""
Manage user authentication with the DivBase server.

This includes login/logout and the getting, storing, using, and refreshing of access + refresh tokens
"""

import time
from dataclasses import dataclass
from pathlib import Path

import httpx
import yaml

# TODO - update user config path here too.
DEFAULT_TOKEN_PATH = Path.home() / ".config" / "divbase" / ".env"


@dataclass
class TokenData:
    """
    Class to hold user token information.
    """

    access_token: str
    refresh_token: str
    access_token_expires_at: int
    refresh_token_expires_at: int

    def dump_tokens(self, output_path: Path) -> None:
        """Dump the user token data to the specified output path"""
        output_path.parent.mkdir(parents=True, exist_ok=True)

        token_dict = {
            "access_token": self.access_token,
            "refresh_token": self.refresh_token,
            "access_token_expires_at": self.access_token_expires_at,
            "refresh_token_expires_at": self.refresh_token_expires_at,
        }
        with open(output_path, "w") as file:
            yaml.safe_dump(token_dict, file, sort_keys=False)

    def is_access_token_expired(self) -> bool:
        """Check if the access token is expired"""
        return time.time() >= self.access_token_expires_at

    def is_refresh_token_expired(self) -> bool:
        """Check if the refresh token is expired"""
        return time.time() >= self.refresh_token_expires_at


def login_to_divbase(email: str, password: str, divbase_url: str) -> None:
    """
    Log in to the DivBase server and return user tokens.
    """
    response = httpx.post(
        f"{divbase_url}/api/v1/auth/login",
        data={
            "grant_type": "password",
            "username": email,  # OAuth2 uses 'username', not 'email'
            "password": password,
        },
        headers={"Content-Type": "application/x-www-form-urlencoded"},
    )
    response.raise_for_status()

    data = response.json()
    token_data = TokenData(
        access_token=data["access_token"],
        refresh_token=data["refresh_token"],
        access_token_expires_at=data["access_token_expires_at"],
        refresh_token_expires_at=data["refresh_token_expires_at"],
    )
    token_data.dump_tokens(output_path=DEFAULT_TOKEN_PATH)


def logout_of_divbase() -> None:
    """
    Log out of the DivBase server.
    TODO - Decide whether to implement token blacklisting on server side.
    """
    if DEFAULT_TOKEN_PATH.exists():
        DEFAULT_TOKEN_PATH.unlink()


def load_user_tokens(token_path: Path = DEFAULT_TOKEN_PATH) -> TokenData:
    """
    Load user tokens from the specified path.
    """
    if not token_path.exists():
        raise FileNotFoundError(f"Your Tokens were not found at {token_path}. Please check you are logged in first.")

    with open(token_path, "r") as file:
        token_dict = yaml.safe_load(file)

    return TokenData(
        access_token=token_dict["access_token"],
        refresh_token=token_dict["refresh_token"],
        access_token_expires_at=token_dict["access_token_expires_at"],
        refresh_token_expires_at=token_dict["refresh_token_expires_at"],
    )
