"""
Authentication for FastAPI routes.

The dependcies here are split in two groups:
1) Deps for routes that use cookies for auth (i.e. the frontend routes)
2) Deps for routes that use direct API access with bearer tokens (i.e. the /api/v1/* routes)
    Auth for API routes can be handled in 2 ways:
    - JWT access token in the Authorization header (for CLI/API clients) logging in via password
    - Personal Access Token (PAT) in the Authorization header (for CLI/API clients) sending requests via PATs.
    PATs cannot be used on the frontend routes.

Several of the dependencies/functions in this file rely on a logged in user, handled by:
- get_current_user for direct API access routes,
- get_current_user_from_cookie for frontend based routes.

The dependencies that depend on this can use e.g. get_current_user as a sub-dependency,
so all checks in the sub-dependency function are ran when you use this dependency.
(see here: https://fastapi.tiangolo.com/yo/advanced/security/oauth2-scopes/#dependency-tree-and-scopes)
"""

import logging
from typing import Annotated

from fastapi import Cookie, Depends, Response
from fastapi.security import OAuth2PasswordBearer
from pydantic import SecretStr, ValidationError
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.auth import verify_user_from_access_token, verify_user_from_refresh_token
from divbase_api.crud.personal_access_tokens import verify_user_from_personal_access_token
from divbase_api.crud.projects import (
    get_project_by_name_or_id_with_user_role,
    pat_role_above_user_role,
)
from divbase_api.db import get_db
from divbase_api.exceptions import AuthenticationError, AuthorizationError
from divbase_api.models.projects import ProjectDB, ProjectRoles
from divbase_api.models.users import UserDB
from divbase_api.security import PAT_TOKEN_PREFIX, TokenType, create_token
from divbase_lib.api_schemas.personal_access_tokens import PATPermissions

logger = logging.getLogger(__name__)

# Returned for JWT tokens and PATs which are not scoped.
# Effectively means the token has whatever permissions the user has.
PAT_NOT_SCOPED = PATPermissions(all_projects=True, task_history=True)

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
        user = await verify_user_from_access_token(db=db, token=access_token)
        if user:
            return user

    # now try refresh token, if this is valid, we will give the user a new access token
    if not refresh_token:
        return None

    user = await verify_user_from_refresh_token(db=db, token=refresh_token)
    if not user:
        return None

    if response:
        new_token_data = create_token(subject=user.id, token_type=TokenType.ACCESS)
        response.set_cookie(
            key=TokenType.ACCESS.value,
            value=new_token_data.token,
            expires=new_token_data.expires_at,
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
    project, user_role = await get_project_by_name_or_id_with_user_role(
        db=db, project_id=project_id, user_id=current_user.id
    )

    if user_role not in ProjectRoles:
        raise AuthorizationError("You don't have permission to access this project.")

    return project, current_user, user_role


#### Dependencies for direct API access routes (using bearer tokens) below ####


async def get_current_user(
    token: Annotated[str, Depends(oauth2_scheme)], db: AsyncSession = Depends(get_db)
) -> tuple[UserDB, PATPermissions]:
    """
    Authenticate from a bearer token (JWT access token or PAT) and return the user plus their scopes.

    The scopes are always a PATPermissions instance:
    - PAT_NOT_SCOPED for JWT tokens (no PAT restrictions apply).
    - The PAT's own PATPermissions for scoped PATs.

    If wondering why no refresh token handling here like in the frontend:
    - The CLI can send the refresh token to "auth/refresh" endpoint instead to get a new access token.
    - PATs do not use refresh tokens.
    """
    # PAT path
    if token.startswith(PAT_TOKEN_PREFIX):
        result = await verify_user_from_personal_access_token(db=db, raw_pat=SecretStr(token))
        if result is None:
            raise AuthenticationError("Invalid or expired personal access token")

        user, pat = result
        try:
            scopes = PATPermissions.model_validate(pat.permissions)
        except ValidationError as err:
            logger.error(
                f"Unexpected error validating PAT permissions for user id={user.id} with permissions dict={pat.permissions}: {err}"
            )
            raise AuthorizationError("Invalid permissions in personal access token.") from None
        return user, scopes

    # JWT path
    user = await verify_user_from_access_token(db=db, token=token)
    if user is None:
        raise AuthenticationError("Invalid or expired access token")
    return user, PAT_NOT_SCOPED


async def get_current_admin_user(
    current_user_and_scopes: Annotated[tuple[UserDB, PATPermissions], Depends(get_current_user)],
) -> UserDB:
    """
    Verify current user has admin access.

    Note: This function works by using a sub-dependency to get the current user,
    so all the checks in the sub-dependency function are ran when you use this dependency.
    (see here: https://fastapi.tiangolo.com/yo/advanced/security/oauth2-scopes/#dependency-tree-and-scopes)
    """
    current_user, _ = current_user_and_scopes
    if not current_user.is_admin:
        raise AuthorizationError("Admin access required")
    return current_user


async def get_project_member(
    project_name: str,
    current_user_and_scopes: Annotated[tuple[UserDB, PATPermissions], Depends(get_current_user)],
    db: AsyncSession = Depends(get_db),
) -> tuple[ProjectDB, UserDB, ProjectRoles]:
    """
    Validates that the user is a member of the project and returns (project, user, effective_role).

    For PATs with scopes, the effective role is the lesser of the user's membership role and
    the role specified in the PAT's project scope. PATs with scopes that don't include the
    requested project will be denied outright.
    """
    current_user, scopes = current_user_and_scopes

    project, user_role = await get_project_by_name_or_id_with_user_role(
        db=db, project_name=project_name, user_id=current_user.id
    )

    if user_role not in ProjectRoles:
        raise AuthorizationError("You don't have permission to access this project.")

    if scopes.all_projects:
        # JWT token or a PAT with no project scoping — use the user's actual role.
        return project, current_user, user_role

    pat_role_str = scopes.projects.get(str(project.id))
    if pat_role_str is None:
        raise AuthorizationError("This personal access token does not have access to this project.")

    try:
        pat_role = ProjectRoles(pat_role_str.lower())
    except ValueError as err:
        logger.error(
            f"Unexpected error validating PAT role for user id={current_user.id} with invalid role='{pat_role_str}' in scopes.projects for project id={project.id}: {err}"
        )
        raise AuthorizationError(f"Invalid role '{pat_role_str}' in personal access token permissions.") from None

    effective_role = pat_role
    if pat_role_above_user_role(pat_role=pat_role, user_role=user_role):
        effective_role = user_role

    return project, current_user, effective_role


async def require_task_history_scope(
    current_user_and_scopes: Annotated[tuple[UserDB, PATPermissions], Depends(get_current_user)],
) -> UserDB:
    """
    Verify the bearer token has the 'task_history' scope.

    JWT tokens (represented by PAT_NOT_SCOPED) always pass.
    Scoped PATs must have task_history=True.
    """
    current_user, scopes = current_user_and_scopes
    if not scopes.task_history:
        raise AuthorizationError("This personal access token does not have the 'task_history' scope.")
    return current_user
