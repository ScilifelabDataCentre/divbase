"""
Authentication endpoints for login and registration.

This is for CLI/API clients to login and get tokens for authenticated requests.
Frontend routes are located in TODO...

TODO:
- Currently only id appended to token payload, could add email + is_admin - but this is debated what is best practice.
"""

import logging

from fastapi import APIRouter, Depends, status
from fastapi.security import OAuth2PasswordRequestForm
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.auth import authenticate_user, user_account_valid
from divbase_api.crud.users import get_user_by_id
from divbase_api.db import get_db
from divbase_api.deps import get_current_user
from divbase_api.exceptions import AuthenticationError
from divbase_api.models.users import UserDB
from divbase_api.schemas.auth import (
    CLILoginResponse,
    RefreshTokenRequest,
    RefreshTokenResponse,
)
from divbase_api.schemas.users import UserResponse
from divbase_api.security import TokenType, create_access_token, create_refresh_token, verify_token

logger = logging.getLogger(__name__)

auth_router = APIRouter()


@auth_router.post("/login", response_model=CLILoginResponse, status_code=status.HTTP_200_OK)
async def login_endpoint(form_data: OAuth2PasswordRequestForm = Depends(), db: AsyncSession = Depends(get_db)):
    """
    Login endpoint to authenticate a user and give them an access + refresh token.

    This is used by CLI/API clients to login and get tokens for authenticated requests.
    Frontend route will has a different endpoint.
    """
    # (form_data.username is defined by OAuth2PasswordRequestForm, but using email for this is fine)
    user = await authenticate_user(db, email=form_data.username, password=form_data.password)
    logger.info(f"User {user.email} logged in successfully.")

    access_token, exp_access = create_access_token(subject=user.id)
    refresh_token, exp_refresh = create_refresh_token(subject=user.id)

    return CLILoginResponse(
        access_token=access_token,
        access_token_expires_at=exp_access,
        refresh_token=refresh_token,
        refresh_token_expires_at=exp_refresh,
        email=user.email,
    )


@auth_router.post("/refresh", response_model=RefreshTokenResponse, status_code=status.HTTP_200_OK)
async def refresh_token_endpoint(refresh_token: RefreshTokenRequest, db: AsyncSession = Depends(get_db)):
    """
    Refresh token endpoint.

    Creates a new access token using the refresh token.
    Unlike the access token pathway, refresh token validates a user account is still valid (active, not deleted, email verified).
    """
    user_id = verify_token(token=refresh_token.refresh_token, desired_token_type=TokenType.REFRESH)
    if not user_id:
        raise AuthenticationError(
            message="Invalid refresh token, please log in again",
        )

    user = await get_user_by_id(db=db, id=user_id)
    if not user:
        logger.warning(f"Attempt to use refresh token for a non-existent user id: {user_id}")
        raise AuthenticationError(message="User not found or inactive or deleted")
    if not user_account_valid(user):
        logger.warning(f"Attempt to refresh token for invalid user account: {user.email} (id: {user.id})")
        raise AuthenticationError(message="User not found or inactive or deleted")

    access_token, expires_at = create_access_token(subject=user.id)
    return RefreshTokenResponse(access_token=access_token, expires_at=expires_at)


@auth_router.post("/logout", status_code=status.HTTP_200_OK)
async def logout_endpoint():
    # TODO
    pass


@auth_router.get("/whoami", status_code=status.HTTP_200_OK, response_model=UserResponse)
async def whoami_endpoint(current_user: UserDB = Depends(get_current_user)):
    """Endpoint to return current logged in user's details."""
    return UserResponse.model_validate(current_user)
