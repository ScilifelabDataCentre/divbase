"""
Authentication endpoints for login and registration.

This is for CLI/API clients to login and get tokens for authenticated requests.
Frontend routes are located in TODO...

TODO:
- Currently only id appended to token payload, could add email + is_admin - but this is debated what is best practice.
"""

import logging

from fastapi import APIRouter, Depends, HTTPException, status
from fastapi.security import OAuth2PasswordRequestForm
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.auth import authenticate_user
from divbase_api.crud.users import get_user_by_id
from divbase_api.db import get_db
from divbase_api.schemas.auth import (
    CLILoginResponse,
    RefreshTokenRequest,
    RefreshTokenResponse,
)
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
    if not user:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid credentials",
            headers={"WWW-Authenticate": "Bearer"},
        )
    if not user.is_active:
        logger.warning(f"Inactive user {user.email} attempted to log in.")
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Inactive user",
        )

    logger.info(f"User {user.email} logged in successfully.")
    return CLILoginResponse(
        access_token=create_access_token(subject=user.id),
        refresh_token=create_refresh_token(subject=user.id),
        email=user.email,
    )


@auth_router.post("/refresh", response_model=RefreshTokenResponse, status_code=status.HTTP_200_OK)
async def refresh_token_endpoint(refresh_token: RefreshTokenRequest, db: AsyncSession = Depends(get_db)):
    """
    Refresh token endpoint.

    Creates a new access token using the refresh token.
    Unlike access token pathway, refresh token validates a user still exists and is still active.

    TODO: Decided if refresh token should not also be refreshed here.
    """

    user_id = verify_token(token=refresh_token.refresh_token, desired_token_type=TokenType.REFRESH)
    if not user_id:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid refresh token, please log in again",
            headers={"WWW-Authenticate": "Bearer"},
        )

    user = await get_user_by_id(db=db, id=user_id)
    if not user or not user.is_active:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="User not found or inactive",
            headers={"WWW-Authenticate": "Bearer"},
        )

    access_token = create_access_token(subject=user.id)
    return RefreshTokenResponse(access_token=access_token)


@auth_router.post("/logout", status_code=status.HTTP_200_OK)
async def logout_endpoint():
    # TODO
    pass
