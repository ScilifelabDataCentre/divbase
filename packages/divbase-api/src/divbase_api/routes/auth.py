"""
Authentication endpoints for login and registration.
"""

from fastapi import APIRouter, Depends, HTTPException, status
from fastapi.security import OAuth2PasswordRequestForm
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.auth import authenticate_user
from divbase_api.db import get_db

auth_router = APIRouter()


@auth_router.post("/login")
async def login_endpoint(form_data: OAuth2PasswordRequestForm = Depends(), db: AsyncSession = Depends(get_db)):
    """
    Login endpoint to authenticate a user and return an access + refresh token (TODO).

    Access + refresh token logic is TODO.
    """
    # (email is username in OAuth2PasswordRequestForm)
    user = await authenticate_user(db, email=form_data.username, password=form_data.password)
    if not user:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid credentials",
            headers={"WWW-Authenticate": "Bearer"},
        )
    return {"message": "Login successful - token logic is TODO"}
