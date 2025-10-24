"""
Security utilities for handling passwords and (TODO) JWT tokens.

bcrypt recommended by FastAPI docs: https://fastapi.tiangolo.com/yo/tutorial/security/oauth2-jwt/
Follows setup from official full stack template: https://github.com/fastapi/full-stack-fastapi-template/blob/master/backend/app/core/security.py
"""

from datetime import datetime, timedelta, timezone
from enum import Enum
from typing import Any

import jwt
from pwdlib import PasswordHash
from pwdlib.hashers.argon2 import Argon2Hasher
from pydantic import SecretStr

from divbase_api.config import settings

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


def create_access_token(subject: str | Any) -> tuple[str, int]:
    """
    Create a new access token for a user. Return token + expiration timestamp.

    Token types differ by the "type" field in the payload, which is either "access" or "refresh".
    """
    expire = datetime.now(timezone.utc) + timedelta(seconds=settings.jwt.access_token_expires_seconds)
    to_encode = {"exp": expire, "sub": str(subject), "type": TokenType.ACCESS.value}
    encoded_jwt = jwt.encode(to_encode, settings.jwt.secret_key.get_secret_value(), algorithm=settings.jwt.algorithm)
    return encoded_jwt, int(expire.timestamp())


def create_refresh_token(subject: str | Any) -> tuple[str, int]:
    """
    Create a new refresh token for a user. Return token + expiration timestamp.

    Token types differ by the "type" field in the payload, which is either "access" or "refresh".
    """
    expire = datetime.now(timezone.utc) + timedelta(seconds=settings.jwt.refresh_token_expires_seconds)
    to_encode = {"exp": expire, "sub": str(subject), "type": TokenType.REFRESH.value}
    encoded_jwt = jwt.encode(to_encode, settings.jwt.secret_key.get_secret_value(), algorithm=settings.jwt.algorithm)
    return encoded_jwt, int(expire.timestamp())


def create_email_verification_token(subject: str | Any) -> tuple[str, int]:
    """
    Create a new email verification token for a user. Return token + expiration timestamp.
    """
    expire = datetime.now(timezone.utc) + timedelta(seconds=settings.email.email_verify_expires_seconds)
    to_encode = {"exp": expire, "sub": str(subject), "type": TokenType.EMAIL_VERIFICATION.value}
    encoded_jwt = jwt.encode(to_encode, settings.jwt.secret_key.get_secret_value(), algorithm=settings.jwt.algorithm)
    return encoded_jwt, int(expire.timestamp())


def create_password_reset_token(subject: str | Any) -> tuple[str, int]:
    """
    Create a new password reset token for a user. Return token + expiration timestamp.
    """
    expire = datetime.now(timezone.utc) + timedelta(seconds=settings.email.password_reset_expires_seconds)
    to_encode = {"exp": expire, "sub": str(subject), "type": TokenType.PASSWORD_RESET.value}
    encoded_jwt = jwt.encode(to_encode, settings.jwt.secret_key.get_secret_value(), algorithm=settings.jwt.algorithm)
    return encoded_jwt, int(expire.timestamp())


def verify_token(token: str, desired_token_type: TokenType) -> int | None:
    """
    Verify and decode JWT token. If successful return the user id, else None.

    If verification fails we should pass not any error information/context back.

    Expiration time is automatically verified in jwt.decode() -> raises jwt.ExpiredSignatureError
    """
    try:
        payload = jwt.decode(
            jwt=token, key=settings.jwt.secret_key.get_secret_value(), algorithms=[settings.jwt.algorithm]
        )
    except (jwt.ExpiredSignatureError, jwt.InvalidTokenError):
        return None

    if payload.get("type") != desired_token_type.value:
        return None
    return int(payload.get("sub"))
