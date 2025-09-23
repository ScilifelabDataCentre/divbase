"""
Security utilities for handling passwords and (TODO) JWT tokens.

bcrypt recommended by FastAPI docs: https://fastapi.tiangolo.com/yo/tutorial/security/oauth2-jwt/
Follows setup from official full stack template: https://github.com/fastapi/full-stack-fastapi-template/blob/master/backend/app/core/security.py
"""

from datetime import datetime, timedelta, timezone
from enum import Enum
from typing import Any

import jwt
from passlib.context import CryptContext
from pydantic import SecretStr

from divbase_api.config import settings

pwd_context = CryptContext(schemes=["bcrypt"], deprecated="auto")


def verify_password(plain_password: str, hashed_password: str) -> bool:
    return pwd_context.verify(plain_password, hashed_password)


def get_password_hash(password: SecretStr) -> str:
    return pwd_context.hash(password.get_secret_value())


class TokenType(str, Enum):
    """Types of JWT tokens used for auth"""

    ACCESS = "access_token"
    REFRESH = "refresh_token"


def create_access_token(subject: str | Any) -> str:
    """
    Create a new access token for a user.

    Token types differ by the "type" field in the payload, which is either "access" or "refresh".
    """
    expire = datetime.now(timezone.utc) + timedelta(seconds=settings.jwt.access_token_expires_seconds)
    to_encode = {"exp": expire, "sub": str(subject), "type": TokenType.ACCESS.value}
    encoded_jwt = jwt.encode(to_encode, settings.jwt.secret_key.get_secret_value(), algorithm=settings.jwt.algorithm)
    return encoded_jwt


def create_refresh_token(subject: str | Any) -> str:
    """
    Create a new refresh token for a user.

    Token types differ by the "type" field in the payload, which is either "access" or "refresh".
    """
    expire = datetime.now(timezone.utc) + timedelta(seconds=settings.jwt.refresh_token_expires_seconds)
    to_encode = {"exp": expire, "sub": str(subject), "type": TokenType.REFRESH.value}
    encoded_jwt = jwt.encode(to_encode, settings.jwt.secret_key.get_secret_value(), algorithm=settings.jwt.algorithm)
    return encoded_jwt


def verify_token(token: str, desired_token_type: TokenType) -> int | None:
    """
    Verify and decode JWT token. If succesful return the user id, else None.

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
