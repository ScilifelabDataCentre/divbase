"""
Authentication for FastAPI routes.

The dependcies here are split in two groups:
1) Deps for routes that use cookies for auth (i.e. the frontend routes)
2) Deps for routes that use direct API access with bearer tokens (i.e. the /api/v1/* routes)

Several of the dependencies/functions in this file rely on a logged in user, handled by:
- get_current_user for direct API access routes,
- get_current_user_from_cookie for frontend based routes.

The dependencies that depend on this can use e.g. get_current_user as a sub-dependency,
so all checks in the sub-dependancy function are ran when you use this dependency.
(see here: https://fastapi.tiangolo.com/yo/advanced/security/oauth2-scopes/#dependency-tree-and-scopes)
"""

import logging
from typing import Annotated

from fastapi import Cookie, Depends, Response
from fastapi.security import OAuth2PasswordBearer
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.auth import user_account_valid
from divbase_api.crud.projects import get_project_id_from_name, get_project_with_user_role
from divbase_api.crud.users import get_user_by_id
from divbase_api.db import get_db
from divbase_api.exceptions import AuthenticationError, AuthorizationError, ProjectNotFoundError
from divbase_api.models.projects import ProjectDB, ProjectRoles
from divbase_api.models.users import UserDB
from divbase_api.security import TokenType, create_token, verify_token

logger = logging.getLogger(__name__)

oauth2_scheme = OAuth2PasswordBearer(tokenUrl="/api/v1/auth/login")


async def _authenticate_frontend_user_from_tokens(
    access_token: str | None,
    refresh_token: str | None,
    db: AsyncSession,
    response: Response | None = None,
) -> UserDB | None:
    """
    Helper function to authenticate a user from either their access token
    (or if that's expired) their refresh token.

    Used for frontend routes where both access + refresh tokens are sent in all requests
    as HttpOnly cookies.
    """
    if not access_token and not refresh_token:
        return None

    if access_token:
        user_id = verify_token(token=access_token, desired_token_type=TokenType.ACCESS)
        if user_id:
            user = await get_user_by_id(db=db, id=user_id)
            if not user:
                return None
            return user

    # now try refresh token, if this is valid, we will give the user a new access token
    if not refresh_token:
        return None

    user_id = verify_token(token=refresh_token, desired_token_type=TokenType.REFRESH)
    if not user_id:
        return None
    user = await get_user_by_id(db=db, id=user_id)
    if not user:
        return None
    if not user_account_valid(user):
        return None

    if response:
        new_access_token, expires_at = create_token(subject=user.id, token_type=TokenType.ACCESS)
        response.set_cookie(
            key=TokenType.ACCESS.value,
            value=new_access_token,
            expires=expires_at,
            httponly=True,
            secure=True,
            samesite="lax",
        )
    return user


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
    return await _authenticate_frontend_user_from_tokens(
        access_token=access_token, refresh_token=refresh_token, db=db, response=response
    )


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
    user = await _authenticate_frontend_user_from_tokens(
        access_token=access_token, refresh_token=refresh_token, db=db, response=response
    )
    if not user:
        raise AuthenticationError("Authentication required")
    return user


async def get_project_member_from_cookie(
    project_id: int,
    current_user: Annotated[UserDB, Depends(get_current_user_from_cookie)],
    db: AsyncSession = Depends(get_db),
) -> tuple[ProjectDB, UserDB, ProjectRoles]:
    """
    This function checks that the user is a member of the project and
    returns the project, user and their role.
    """
    project, user_role = await get_project_with_user_role(db=db, project_id=project_id, user_id=current_user.id)

    if not project:
        raise ProjectNotFoundError()

    if user_role not in ProjectRoles:
        raise AuthorizationError("You don't have permission to access this project.")

    return project, current_user, user_role


async def get_current_user(token: Annotated[str, Depends(oauth2_scheme)], db: AsyncSession = Depends(get_db)) -> UserDB:
    """
    Get current user from JWT Access token

    We should not be specific about why/if credentials are invalid.

    If wondering why no refresh token handling here like in the frontend:
    - The CLI can send the refresh token to "auth/refresh" endpoint instead to get a new access token.
    """
    user_id = verify_token(token, desired_token_type=TokenType.ACCESS)
    if user_id is None:
        raise AuthenticationError("Authentication required")

    user = await get_user_by_id(db=db, id=user_id)
    if not user or not user.is_active or user.is_deleted:
        raise AuthenticationError("Account does not exist or is inactive or is deleted")
    if not user.email_verified:
        raise AuthenticationError(
            "Email address not verified, check your inbox or visit the website to resend a verification email."
        )

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


async def get_project_member(
    project_name: str,
    current_user: Annotated[UserDB, Depends(get_current_user)],
    db: AsyncSession = Depends(get_db),
) -> tuple[ProjectDB, UserDB, ProjectRoles]:
    """
    This function checks that the user is a member of the project and
    returns the project, user and their role.
    """
    project_id = await get_project_id_from_name(db=db, project_name=project_name)
    if not project_id:
        raise ProjectNotFoundError()

    project, user_role = await get_project_with_user_role(db=db, project_id=project_id, user_id=current_user.id)
    if not project:
        raise ProjectNotFoundError()

    if user_role not in ProjectRoles:
        raise AuthorizationError("You don't have permission to access this project.")

    return project, current_user, user_role
