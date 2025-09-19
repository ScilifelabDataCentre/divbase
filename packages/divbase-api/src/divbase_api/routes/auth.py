"""
Authentication endpoints for login and registration.

This is for CLI/API clients to login and get tokens for authenticated requests.
Frontend routes are located in TODO...

TODO:
- Currently only id appended to token payload, could add email + is_admin - but this is debated what is best practice.
"""

from fastapi import APIRouter, Depends, HTTPException, status
from fastapi.security import OAuth2PasswordRequestForm
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.auth import authenticate_user
from divbase_api.db import get_db
from divbase_api.schemas.auth import CLILoginResponse
from divbase_api.security import create_access_token, create_refresh_token

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
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Inactive user",
        )

    return CLILoginResponse(
        access_token=create_access_token(subject=user.id),
        refresh_token=create_refresh_token(subject=user.id),
        email=user.email,
    )
