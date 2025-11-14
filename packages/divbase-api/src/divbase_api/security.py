"""
Security utilities for handling passwords and JSON Web Tokens (JWTs).

FastAPI recommend using pwdlib and argon2 for password hashing (https://fastapi.tiangolo.com/tutorial/security/oauth2-jwt/)

Follows setup from official full stack template: https://github.com/fastapi/full-stack-fastapi-template/blob/master/backend/app/core/security.py
"""

from datetime import datetime, timedelta, timezone
from enum import Enum
from typing import Any

import jwt
from pwdlib import PasswordHash
from pwdlib.hashers.argon2 import Argon2Hasher
from pydantic import SecretStr

from divbase_api.api_config import settings

password_hash = PasswordHash(hashers=[Argon2Hasher()])


def verify_password(plain_password: str, hashed_password: str) -> bool:
    return password_hash.verify(plain_password, hashed_password)


def get_password_hash(password: SecretStr) -> str:
    return password_hash.hash(password.get_secret_value())


class TokenType(str, Enum):
    """Types of JWT tokens used for e.g. auth or email verification/password reset."""

    ACCESS = "access_token"
    REFRESH = "refresh_token"
    EMAIL_VERIFICATION = "email_verification"
    PASSWORD_RESET = "password_reset"


token_expires_delta: dict[TokenType, timedelta] = {
    TokenType.ACCESS: timedelta(seconds=settings.jwt.access_token_expires_seconds),
    TokenType.REFRESH: timedelta(seconds=settings.jwt.refresh_token_expires_seconds),
    TokenType.EMAIL_VERIFICATION: timedelta(seconds=settings.email.email_verify_expires_seconds),
    TokenType.PASSWORD_RESET: timedelta(seconds=settings.email.password_reset_expires_seconds),
}


def create_token(subject: str | Any, token_type: TokenType) -> tuple[str, int]:
    """
    Create a JWT token for access, refresh, email verification or password reset.

    Returns a tuple of the JWT and expiry UNIX time stamp for the token.

    Tokens specify the "type" field in the payload,
    so we can validate an access token is not used for e.g. password reset.
    """
    expire = datetime.now(timezone.utc) + token_expires_delta[token_type]
    to_encode = {"exp": expire, "iat": datetime.now(timezone.utc), "sub": str(subject), "type": token_type.value}
    encoded_jwt = jwt.encode(to_encode, settings.jwt.secret_key.get_secret_value(), algorithm=settings.jwt.algorithm)
    return encoded_jwt, int(expire.timestamp())


def verify_token(token: str, desired_token_type: TokenType) -> int | None:
    """
    Verify and decode JWT token. If successful return the user id, else None.

    If verification fails we should pass not any error information/context back.

    Expiration time is automatically verified in jwt.decode() -> raises jwt.ExpiredSignatureError
    """
    if desired_token_type == TokenType.REFRESH:
        raise ValueError("Use verify_refresh_token() for refresh tokens.")
    try:
        payload = jwt.decode(
            jwt=token, key=settings.jwt.secret_key.get_secret_value(), algorithms=[settings.jwt.algorithm]
        )
    except (jwt.ExpiredSignatureError, jwt.InvalidTokenError):
        return None

    if payload.get("type") != desired_token_type:
        return None
    return int(payload.get("sub"))


def verify_refresh_token(token: str) -> tuple[int, datetime] | None:
    """
    Verify and decode a JWT token of TokenType type "refresh".

    This is implemented separately to normal verify_token() because we also want to return the "iat" (issued at) time.
    For a password_last_updated comparison check.

    If successful return the user id and iat, else return None.
    """
    try:
        payload = jwt.decode(
            jwt=token, key=settings.jwt.secret_key.get_secret_value(), algorithms=[settings.jwt.algorithm]
        )
    except (jwt.ExpiredSignatureError, jwt.InvalidTokenError):
        return None

    if payload.get("type") != TokenType.REFRESH:
        return None

    iat = datetime.fromtimestamp(payload.get("iat"), tz=timezone.utc)
    return int(payload.get("sub")), iat


def verify_expired_token(token: str, desired_token_type: TokenType) -> int | None:
    """
    Verify and decode a (potentially) expired JWT token. If successful return the user id, else None.

    Used e.g. for password reset/email verification flow where we want to inform user that their token has expired,
    or email is already verified. So for UX reasons.
    (We're still checking the signature and token type, just not expiration time.)
    """
    if desired_token_type not in (TokenType.EMAIL_VERIFICATION, TokenType.PASSWORD_RESET):
        raise ValueError("Can only verify expired tokens for email verification or password reset.")

    try:
        payload = jwt.decode(
            jwt=token,
            key=settings.jwt.secret_key.get_secret_value(),
            algorithms=[settings.jwt.algorithm],
            options={"verify_exp": False},
        )
    except jwt.InvalidTokenError:
        return None

    if payload.get("type") != desired_token_type:
        return None
    return int(payload.get("sub"))
