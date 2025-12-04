"""
Security utilities for handling passwords and JSON Web Tokens (JWTs).

FastAPI recommend using pwdlib and argon2 for password hashing (https://fastapi.tiangolo.com/tutorial/security/oauth2-jwt/)

Follows setup from official full stack template: https://github.com/fastapi/full-stack-fastapi-template/blob/master/backend/app/core/security.py
"""

import uuid
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from enum import StrEnum
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


class TokenType(StrEnum):
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


@dataclass
class TokenData:
    """Returned when creating a JWT."""

    token: str
    expires_at: int


@dataclass
class VerifiedTokenData:
    """Returned after verifying and decoding a JWT."""

    user_id: int
    issued_at: datetime
    jti: str


def create_token(subject: str | Any, token_type: TokenType) -> TokenData:
    """
    Create a JWT token for access, refresh, email verification or password reset.

    Returns a tuple of the JWT and expiry UNIX time stamp for the token.

    Tokens specify the "type" field in the payload,
    so we can validate an access token is not used for e.g. password reset.
    """
    expire = datetime.now(timezone.utc) + token_expires_delta[token_type]
    jti = str(uuid.uuid4())

    to_encode = {
        "exp": expire,
        "iat": datetime.now(timezone.utc),
        "jti": jti,
        "sub": str(subject),
        "type": token_type.value,
    }
    encoded_jwt = jwt.encode(to_encode, settings.jwt.secret_key.get_secret_value(), algorithm=settings.jwt.algorithm)
    return TokenData(token=encoded_jwt, expires_at=int(expire.timestamp()))


def verify_token(token: str, desired_token_type: TokenType) -> VerifiedTokenData | None:
    """
    Verify and decode JWT token. If successful, return the verified token data; else None.

    If verification fails we should not pass any error information/context back.

    Expiration time is automatically verified in jwt.decode() -> raises jwt.ExpiredSignatureError
    """
    try:
        payload = jwt.decode(
            jwt=token, key=settings.jwt.secret_key.get_secret_value(), algorithms=[settings.jwt.algorithm]
        )
    except (jwt.ExpiredSignatureError, jwt.InvalidTokenError):
        return None

    if payload.get("type") != desired_token_type:
        return None

    # edge case: if we update the JWT structure and these fields are missing
    # as token issued before the change
    if not all(payload.get(field) for field in ["sub", "iat", "jti"]):
        return None

    return VerifiedTokenData(
        user_id=int(payload["sub"]),
        issued_at=datetime.fromtimestamp(payload["iat"], tz=timezone.utc),
        jti=payload["jti"],
    )
