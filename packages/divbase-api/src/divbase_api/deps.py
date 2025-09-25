"""
Authentication for FastAPI routes.

TODO - perhaps token types should be a str enum
"""

import logging
from typing import Annotated

from fastapi import Cookie, Depends, Response
from fastapi.security import OAuth2PasswordBearer
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.config import settings
from divbase_api.crud.users import get_user_by_id
from divbase_api.db import get_db
from divbase_api.exceptions import AuthenticationError, AuthorizationError
from divbase_api.models.users import UserDB
from divbase_api.security import TokenType, create_access_token, verify_token

logger = logging.getLogger(__name__)

oauth2_scheme = OAuth2PasswordBearer(tokenUrl="/api/v1/auth/login")


async def get_current_user_from_cookie_optional(
    access_token: str | None = Cookie(None),
    refresh_token: str | None = Cookie(None),
    db: AsyncSession = Depends(get_db),
    response: Response = None,
) -> UserDB | None:
    """
    Get user from the JWT access token stored inside the httponly cookie.

    Returns None if not logged in.
    Use in routes where both logged in and not logged in users are allowed.
    """
    if not access_token and not refresh_token:
        return None

    if access_token:
        user_id = verify_token(token=access_token, desired_token_type=TokenType.ACCESS)
        if user_id:
            user = await get_user_by_id(db=db, id=user_id)
            if not user or not user.is_active:
                return None
            return user

    # now try refresh token, if this is valid, we give user a new access token
    if not refresh_token:
        return None

    user_id = verify_token(token=refresh_token, desired_token_type=TokenType.REFRESH)
    if not user_id:
        return None
    user = await get_user_by_id(db=db, id=user_id)
    if not user or not user.is_active:
        return None

    if response:
        new_access_token = create_access_token(subject=user.id)
        response.set_cookie(
            key=TokenType.ACCESS.value,
            value=new_access_token,
            max_age=settings.jwt.access_token_expires_seconds,
            httponly=True,
            secure=True,
            samesite="lax",
        )
    return user


async def get_current_user_from_cookie(
    access_token: str | None = Cookie(None),
    refresh_token: str | None = Cookie(None),
    db: AsyncSession = Depends(get_db),
    response: Response = None,
) -> UserDB:
    """
    Get user from the JWT access token stored inside the httponly cookie.
    Raises AuthenticationError if not logged in.

    Use in routes where user must be logged in.
    """
    user = await get_current_user_from_cookie_optional(access_token, refresh_token, db, response)
    if not user:
        raise AuthenticationError("Authentication required")
    return user


async def get_current_admin_user_from_cookie(
    current_user: Annotated[UserDB, Depends(get_current_user_from_cookie)],
) -> UserDB:
    """
    Verify current user has admin access.

    Note: This function works by using a sub-dependency to get the current user,
    so all the checks in the sub-dependancy function are ran when you use this dependency.
    (see here: https://fastapi.tiangolo.com/yo/advanced/security/oauth2-scopes/#dependency-tree-and-scopes)
    """
    if not current_user.is_admin:
        raise AuthorizationError("Admin access required")
    return current_user


async def get_current_user(token: Annotated[str, Depends(oauth2_scheme)], db: AsyncSession = Depends(get_db)) -> UserDB:
    """
    Get current user from JWT Access token

    We should not be specific about why/if credentials are invalid.
    """
    user_id = verify_token(token, desired_token_type=TokenType.ACCESS)
    if user_id is None:
        raise AuthenticationError("Authentication required")

    user = await get_user_by_id(db=db, id=user_id)
    if not user or not user.is_active:
        raise AuthenticationError("Account does not exist or is inactive")

    return user


async def get_current_admin_user(current_user: Annotated[UserDB, Depends(get_current_user)]) -> UserDB:
    """
    Verify current user has admin access.

    Note: This function works by using a sub-dependency to get the current user,
    so all the checks in the sub-dependancy function are ran when you use this dependency.
    (see here: https://fastapi.tiangolo.com/yo/advanced/security/oauth2-scopes/#dependency-tree-and-scopes)
    """
    if not current_user.is_admin:
        raise AuthorizationError("Admin access required")
    return current_user
