"""
Frontend routes for authentication-related pages.
"""

from fastapi import APIRouter, Depends, Form, Request
from fastapi.responses import HTMLResponse
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.config import settings
from divbase_api.crud.auth import authenticate_user
from divbase_api.db import get_db
from divbase_api.deps import get_current_user_from_cookie
from divbase_api.frontend_routes.core import templates
from divbase_api.models.users import UserDB
from divbase_api.schemas.users import UserResponse
from divbase_api.security import TokenType, create_access_token, create_refresh_token

fr_auth_router = APIRouter()


@fr_auth_router.get("/login", response_class=HTMLResponse)
async def get_login(request: Request, current_user: UserDB | None = Depends(get_current_user_from_cookie)):
    """Render the login page."""
    if current_user:
        return templates.TemplateResponse(
            request=request,
            name="index.html",
            context={"request": request, "user": UserResponse.model_validate(current_user)},
        )
    return templates.TemplateResponse(request=request, name="login.html")


@fr_auth_router.post("/login", response_class=HTMLResponse)
async def post_login(
    request: Request, email: str = Form(...), password: str = Form(...), db: AsyncSession = Depends(get_db)
):
    """Handle login form submission."""
    user = await authenticate_user(db, email=email, password=password)
    if not user:
        return templates.TemplateResponse(
            request=request,
            name="login.html",
            context={"request": request, "error": "Invalid email or password"},
        )

    access_token = create_access_token(subject=user.id)
    refresh_token = create_refresh_token(subject=user.id)

    response = templates.TemplateResponse(
        request=request,
        name="index.html",
        context={"user": UserResponse.model_validate(user)},
    )

    response.set_cookie(
        key=TokenType.ACCESS.value,
        value=access_token,
        max_age=settings.jwt.access_token_expires_seconds,
        httponly=True,
        secure=False,  # TODO, set to True in Prod.
        samesite="lax",
    )

    response.set_cookie(
        key=TokenType.REFRESH.value,
        value=refresh_token,
        max_age=settings.jwt.refresh_token_expires_seconds,
        httponly=True,
        secure=False,  # TODO, set to True in Prod.
        samesite="lax",
    )

    return response


@fr_auth_router.post("/logout", response_class=HTMLResponse)
async def post_logout(request: Request):
    """Handle logout form submission."""
    # TODO - decide if should store invalid token in DB.
    response = templates.TemplateResponse(
        request=request,
        name="index.html",
        context={"user": None},
    )

    response.delete_cookie(TokenType.ACCESS.value)
    response.delete_cookie(TokenType.REFRESH.value)

    return response
