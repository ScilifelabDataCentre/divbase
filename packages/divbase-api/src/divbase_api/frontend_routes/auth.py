"""
Frontend routes for authentication-related pages.
"""

import logging

from fastapi import APIRouter, BackgroundTasks, Depends, Form, Query, Request, status
from fastapi.responses import HTMLResponse, RedirectResponse
from pydantic import SecretStr
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.api_config import settings
from divbase_api.crud.auth import (
    authenticate_user,
    check_user_email_verified,
    confirm_user_email,
)
from divbase_api.crud.users import create_user, get_user_by_email
from divbase_api.db import get_db
from divbase_api.deps import get_current_user_from_cookie_optional
from divbase_api.exceptions import AuthenticationError
from divbase_api.frontend_routes.core import templates
from divbase_api.models.users import UserDB
from divbase_api.schemas.users import UserCreate, UserResponse
from divbase_api.security import (
    TokenType,
    create_access_token,
    create_refresh_token,
    verify_expired_token,
    verify_token,
)
from divbase_api.services.email_sender import send_verification_email

logger = logging.getLogger(__name__)


fr_auth_router = APIRouter()


@fr_auth_router.get("/login", response_class=HTMLResponse)
async def get_login(request: Request, current_user: UserDB | None = Depends(get_current_user_from_cookie_optional)):
    """
    Render the login page.
    If user is already logged in, redirect to home page.
    """
    if current_user:
        return RedirectResponse(url="/", status_code=status.HTTP_302_FOUND)

    return templates.TemplateResponse(request=request, name="auth_pages/login.html")


@fr_auth_router.post("/login", response_class=HTMLResponse)
async def post_login(
    request: Request, email: str = Form(...), password: str = Form(...), db: AsyncSession = Depends(get_db)
):
    """Handle login form submission."""
    try:
        user = await authenticate_user(db, email=email, password=password)
    except AuthenticationError as e:
        logger.info(f"Failed login attempt for email: {email} - {e.message}")

        return templates.TemplateResponse(
            request=request,
            name="auth_pages/login.html",
            context={"request": request, "error": e.message},
        )

    logger.info(f"User {user.email} logged in successfully via frontend.")
    access_token, access_expires_at = create_access_token(subject=user.id)
    refresh_token, refresh_expires_at = create_refresh_token(subject=user.id)

    response = templates.TemplateResponse(
        request=request,
        name="index.html",
        context={"current_user": UserResponse.model_validate(user)},
    )

    response.set_cookie(
        key=TokenType.ACCESS.value,
        value=access_token,
        expires=access_expires_at,
        httponly=True,
        secure=True,
        samesite="lax",
    )

    response.set_cookie(
        key=TokenType.REFRESH.value,
        value=refresh_token,
        expires=refresh_expires_at,
        httponly=True,
        secure=True,
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
        context={"current_user": None},
    )

    response.delete_cookie(TokenType.ACCESS.value)
    response.delete_cookie(TokenType.REFRESH.value)

    return response


@fr_auth_router.get("/register", response_class=HTMLResponse)
async def get_register(request: Request, current_user: UserDB | None = Depends(get_current_user_from_cookie_optional)):
    """Render the registration page."""
    if current_user:
        return templates.TemplateResponse(
            request=request,
            name="index.html",
            context={"request": request, "current_user": UserResponse.model_validate(current_user)},
        )
    return templates.TemplateResponse(request=request, name="auth_pages/register.html")


@fr_auth_router.post("/register", response_class=HTMLResponse)
async def post_register(
    request: Request,
    background_tasks: BackgroundTasks,
    name: str = Form(...),
    email: str = Form(...),
    password: str = Form(...),
    confirm_password: str = Form(...),
    db: AsyncSession = Depends(get_db),
):
    """Handle registration form submission."""

    def registration_failed_response(error_message: str):
        """Helper to return registration failed response with custom error message."""
        return templates.TemplateResponse(
            request=request,
            name="auth_pages/register.html",
            context={
                "request": request,
                "error": error_message,
                "name": name,
                "email": email,
            },
        )

    if password != confirm_password:
        return registration_failed_response("Passwords do not match")

    existing_user = await get_user_by_email(db=db, email=email)
    if existing_user:  # Not recommended to specify why failed, just say failed.
        return registration_failed_response("Registration failed, please try again.")

    try:
        user_data = UserCreate(name=name, email=email, password=SecretStr(password))
        user = await create_user(db=db, user_data=user_data, is_admin=False)
    except Exception as e:
        logger.error(f"Error creating user: {e}")
        return registration_failed_response("Registration failed, please try again.")

    background_tasks.add_task(send_verification_email, email_to=user.email, user_id=user.id)

    logger.info(f"New user registered: {user_data.email=}")
    return templates.TemplateResponse(
        request=request,
        name="auth_pages/register_success.html",
        context={
            "request": request,
            "name": user_data.name,
            "email": user_data.email,
            "from_email": settings.email.from_email,
        },
    )


@fr_auth_router.get("/verify-email", response_class=HTMLResponse)
async def get_verify_email(
    request: Request,
    token: str = Query(...),
    db: AsyncSession = Depends(get_db),
):
    """
    Handle email verification.

    To access this endpoint a user recieves an email with link to verify their email.
    The link contains a JWT as query param in the URL.

    For the sake of UX, we render different HTML pages depending on the outcome of the verification:
    if token invalid

    """
    expired_token = False

    user_id = verify_token(token=token, desired_token_type=TokenType.EMAIL_VERIFICATION)
    if not user_id:
        expired_token = True
        user_id = verify_expired_token(token=token, desired_token_type=TokenType.EMAIL_VERIFICATION)

    # token is neither valid nor expired
    if not user_id:
        return templates.TemplateResponse(
            request=request,
            name="auth_pages/email_verification.html",
            context={"error": "Invalid email verification link. Please request a new link below."},
        )

    already_verified = await check_user_email_verified(db=db, id=user_id)

    if expired_token and not already_verified:
        return templates.TemplateResponse(
            request=request,
            name="auth_pages/email_verification.html",
            context={"error": "Email verification link has expired. Please request a new link below."},
        )

    if already_verified:
        return templates.TemplateResponse(
            request=request,
            name="auth_pages/login.html",
            context={"info": "Your email has already been verified, you can log in to DivBase directly"},
        )

    # token is valid and email not verified yet proceed to verify email
    user = await confirm_user_email(db=db, id=user_id)
    return templates.TemplateResponse(
        request=request,
        name="auth_pages/login.html",
        context={"success": f"Thank you for verifying your email '{user.email}', you can now log in to DivBase"},
    )


@fr_auth_router.post("/resend-email-verification", response_class=HTMLResponse)
async def resend_verification_email(
    request: Request,
    background_tasks: BackgroundTasks,
    email: str = Form(...),
    db: AsyncSession = Depends(get_db),
):
    """
    Handle resending the email verification link.
    """
    LINK_SENT_MSG = "If your account exists, a verification email has been sent. Please check your inbox."

    user = await get_user_by_email(db=db, email=email)
    if not user:
        # Do not differentiate between existing and non-existing users for security reasons
        return templates.TemplateResponse(
            request=request,
            name="auth_pages/email_verification.html",
            context={"email": email, "success": LINK_SENT_MSG},
        )

    if user.email_verified:
        return templates.TemplateResponse(
            request=request,
            name="auth_pages/login.html",
            context={"info": "Your email has already been verified, you can log in to DivBase directly"},
        )

    background_tasks.add_task(send_verification_email, email_to=user.email, user_id=user.id)

    return templates.TemplateResponse(
        request=request,
        name="auth_pages/email_verification.html",
        context={"email": email, "success": LINK_SENT_MSG},
    )


@fr_auth_router.get("/resend-verification-email", response_class=HTMLResponse)
async def get_resend_verification_email(
    request: Request,
    db: AsyncSession = Depends(get_db),
):
    """
    Display the resend verification email page.
    """
    return templates.TemplateResponse(request=request, name="auth_pages/email_verification.html")
