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

from divbase_api.crud.auth import authenticate_user, verify_user_from_refresh_token
from divbase_api.crud.revoked_tokens import revoke_token_on_logout
from divbase_api.db import get_db
from divbase_api.deps import get_current_user
from divbase_api.exceptions import AuthenticationError
from divbase_api.models.users import UserDB
from divbase_api.schemas.auth import (
    CLILoginResponse,
    LogoutRequest,
    RefreshTokenRequest,
    RefreshTokenResponse,
)
from divbase_api.schemas.users import UserResponse
from divbase_api.security import TokenType, create_token, verify_token

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

    access_token_data = create_token(subject=user.id, token_type=TokenType.ACCESS)
    refresh_token_data = create_token(subject=user.id, token_type=TokenType.REFRESH)

    return CLILoginResponse(
        access_token=access_token_data.token,
        access_token_expires_at=access_token_data.expires_at,
        refresh_token=refresh_token_data.token,
        refresh_token_expires_at=refresh_token_data.expires_at,
        email=user.email,
    )


@auth_router.post("/refresh", response_model=RefreshTokenResponse, status_code=status.HTTP_200_OK)
async def refresh_token_endpoint(refresh_token: RefreshTokenRequest, db: AsyncSession = Depends(get_db)):
    """
    Refresh token endpoint.

    Creates a new access token using the refresh token.
    Unlike the access token pathway, refresh token runs extra validation on a user account:
    (active, not deleted, email verified and password not changed since refresh token issued.)
    """
    user = await verify_user_from_refresh_token(db=db, token=refresh_token.refresh_token)
    if not user:
        raise AuthenticationError(message="Invalid or expired refresh token, please log in again")

    token_data = create_token(subject=user.id, token_type=TokenType.ACCESS)
    return RefreshTokenResponse(access_token=token_data.token, expires_at=token_data.expires_at)


@auth_router.post("/logout", status_code=status.HTTP_200_OK)
async def logout_endpoint(logout_request: LogoutRequest, db: AsyncSession = Depends(get_db)):
    token_data = verify_token(logout_request.refresh_token, TokenType.REFRESH)
    if token_data:
        await revoke_token_on_logout(db=db, token_jti=token_data.jti, user_id=token_data.user_id)


@auth_router.get("/whoami", status_code=status.HTTP_200_OK, response_model=UserResponse)
async def whoami_endpoint(current_user: UserDB = Depends(get_current_user)):
    """Endpoint to return current logged in user's details."""
    return UserResponse.model_validate(current_user)
