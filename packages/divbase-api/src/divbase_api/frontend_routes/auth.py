"""
Frontend routes for authentication-related pages.
"""

import logging

from fastapi import APIRouter, BackgroundTasks, Depends, Form, Query, Request, status
from fastapi.responses import HTMLResponse, RedirectResponse
from pydantic import SecretStr, ValidationError
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.api_config import settings
from divbase_api.crud.auth import (
    authenticate_user,
    check_user_email_verified,
    confirm_user_email,
    delete_auth_cookies,
    update_user_password,
)
from divbase_api.crud.users import create_user, get_user_by_email, get_user_by_id_or_raise
from divbase_api.db import get_db
from divbase_api.deps import get_current_user_from_cookie_optional
from divbase_api.exceptions import AuthenticationError
from divbase_api.frontend_routes.core import templates
from divbase_api.models.users import UserDB
from divbase_api.schemas.users import UserCreate, UserPasswordUpdate, UserResponse
from divbase_api.security import (
    TokenType,
    create_token,
    verify_token,
)
from divbase_api.services.email_sender import (
    send_email_already_verified_email,
    send_password_has_been_reset_email,
    send_password_reset_email,
    send_verification_email,
)

logger = logging.getLogger(__name__)


fr_auth_router = APIRouter()


INVALID_EXPIRED_PASSWORD_TOKEN_MSG = "Invalid or expired reset password link. Please request a new link below."
INVALID_EXPIRED_EMAIL_TOKEN_MSG = "Invalid or expired email verification link. Please request a new link below."


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
    access_token, access_expires_at = create_token(subject=user.id, token_type=TokenType.ACCESS)
    refresh_token, refresh_expires_at = create_token(subject=user.id, token_type=TokenType.REFRESH)

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
    response = RedirectResponse(url="/", status_code=status.HTTP_302_FOUND)
    return delete_auth_cookies(response=response)


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
    current_user: UserDB | None = Depends(get_current_user_from_cookie_optional),
):
    """
    Handle email verification and redirect to confirmation page.

    To access this endpoint a user receives an email with link to verify their email.
    The link contains a JWT as query param in the URL.

    We don't verify email on clicking the link as some mailboxes may click all links automatically for security scanning.
    Instead, we show a confirmation page where a user has to click a button to confirm verification of their email.
    """
    if current_user:
        return RedirectResponse(url="/", status_code=status.HTTP_302_FOUND)

    user_id = verify_token(token=token, desired_token_type=TokenType.EMAIL_VERIFICATION)
    if not user_id:
        return templates.TemplateResponse(
            request=request,
            name="auth_pages/email_verification.html",
            context={"error": INVALID_EXPIRED_EMAIL_TOKEN_MSG},
        )

    user = await get_user_by_id_or_raise(db=db, id=user_id)
    return templates.TemplateResponse(
        request=request,
        name="auth_pages/email_verification_confirm.html",
        context={"token": token, "email": user.email},
    )


@fr_auth_router.post("/confirm-email-verification", response_class=HTMLResponse)
async def confirm_email_verification(
    request: Request,
    token: str = Form(...),
    db: AsyncSession = Depends(get_db),
):
    """
    Confirm email verification after user explicitly clicks a button.
    """
    user_id = verify_token(token=token, desired_token_type=TokenType.EMAIL_VERIFICATION)
    if not user_id:
        return templates.TemplateResponse(
            request=request,
            name="auth_pages/email_verification.html",
            context={"error": INVALID_EXPIRED_EMAIL_TOKEN_MSG},
        )

    already_verified = await check_user_email_verified(db=db, id=user_id)
    if already_verified:
        return templates.TemplateResponse(
            request=request,
            name="auth_pages/login.html",
            context={"info": "Your email has already been verified, you can log in to DivBase directly."},
        )

    user = await confirm_user_email(db=db, id=user_id)
    return templates.TemplateResponse(
        request=request,
        name="auth_pages/login.html",
        context={"success": f"Thank you for verifying your email '{user.email}', you can now log in to DivBase."},
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
        # User is already verified, inform them by email
        # To prevent information leakage (which accounts exist and don't exists),
        # we show the same success message on the frontend, but email them to inform them they can already login.

        background_tasks.add_task(send_email_already_verified_email, email_to=user.email)
        return templates.TemplateResponse(
            request=request,
            name="auth_pages/login.html",
            context={"success": LINK_SENT_MSG},
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
    current_user: UserDB | None = Depends(get_current_user_from_cookie_optional),
):
    """
    Display the resend verification email page.
    """
    # redirect to home if already logged in.
    if current_user:
        return RedirectResponse(url="/", status_code=status.HTTP_302_FOUND)

    return templates.TemplateResponse(request=request, name="auth_pages/email_verification.html")


@fr_auth_router.get("/forgot-password", response_class=HTMLResponse)
async def get_forgot_password_page(
    request: Request,
    db: AsyncSession = Depends(get_db),
    current_user: UserDB | None = Depends(get_current_user_from_cookie_optional),
):
    """
    Display the forgot password page.
    """
    # TODO - think about how to handle logged in user? - Log them out?
    if current_user:
        return RedirectResponse(url="/", status_code=status.HTTP_302_FOUND)

    return templates.TemplateResponse(request=request, name="auth_pages/forgot_password.html")


@fr_auth_router.post("/forgot-password", response_class=HTMLResponse)
async def post_forgot_password_form(
    request: Request,
    background_tasks: BackgroundTasks,
    email: str = Form(...),
    db: AsyncSession = Depends(get_db),
    current_user: UserDB | None = Depends(get_current_user_from_cookie_optional),
):
    """Handle forgot password form submission to send a password reset email."""
    if current_user:
        # TODO - think about how to handle logged in user? - Log them out?
        return RedirectResponse(url="/", status_code=status.HTTP_302_FOUND)

    RESET_LINK_SENT_MSG = (
        f"If your account exists, and your email is verified, a password reset email has been sent to {email}. \n Please check your inbox."
        + f"\nThe email will be sent from {settings.email.from_email}."
    )

    user = await get_user_by_email(db=db, email=email)

    # do not differentiate between existing and non-existing users for security reasons
    if not user or not user.email_verified:
        logger.info(
            f"A password reset email was requested for '{email}' but not sent. User exists: {bool(user)}, email verified: {user.email_verified if user else 'N/A'}"
        )
        return templates.TemplateResponse(
            request=request,
            name="auth_pages/login.html",
            context={"success": RESET_LINK_SENT_MSG},
        )

    background_tasks.add_task(send_password_reset_email, email_to=user.email, user_id=user.id)

    logger.info(f"Password reset email sent to: {email}")
    return templates.TemplateResponse(
        request=request,
        name="auth_pages/login.html",
        context={"success": RESET_LINK_SENT_MSG},
    )


@fr_auth_router.get("/reset-password", response_class=HTMLResponse)
async def get_reset_password_page(
    request: Request,
    token: str,
    db: AsyncSession = Depends(get_db),
    current_user: UserDB | None = Depends(get_current_user_from_cookie_optional),
):
    """
    Display the reset password page.

    To access this endpoint a user receives an email with link to reset their password.
    The link contains a JWT as query param in the URL.
    """
    # TODO - add email to context to show which email is being reset?
    user_id = verify_token(token=token, desired_token_type=TokenType.PASSWORD_RESET)
    if not user_id:
        return templates.TemplateResponse(
            request=request,
            name="auth_pages/forgot_password.html",
            context={"error": INVALID_EXPIRED_PASSWORD_TOKEN_MSG},
        )

    user = await get_user_by_id_or_raise(db=db, id=user_id)
    return templates.TemplateResponse(
        request=request,
        name="auth_pages/reset_password.html",
        context={"token": token, "email": user.email},
    )


@fr_auth_router.post("/reset-password", response_class=HTMLResponse)
async def post_reset_password_form(
    request: Request,
    background_tasks: BackgroundTasks,
    token: str = Form(...),
    password: str = Form(...),
    confirm_password: str = Form(...),
    db: AsyncSession = Depends(get_db),
):
    """
    Handle reset password form submission.
    """
    user_id = verify_token(token=token, desired_token_type=TokenType.PASSWORD_RESET)
    if not user_id:
        return templates.TemplateResponse(
            request=request,
            name="auth_pages/forgot_password.html",
            context={"error": INVALID_EXPIRED_PASSWORD_TOKEN_MSG},
        )

    # Client side validation should mean these are never raised, but always have to check on server side.
    try:
        password_data = UserPasswordUpdate(password=SecretStr(password), confirm_password=SecretStr(confirm_password))
    except ValidationError as err:
        error_msg = err.errors()[0]["msg"]
        if "Value error, " in error_msg:
            error_msg = error_msg.replace("Value error, ", "")

        return templates.TemplateResponse(
            request=request,
            name="auth_pages/reset_password.html",
            context={"request": request, "token": token, "error": str(error_msg)},
        )
    user = await update_user_password(db=db, user_id=user_id, password_data=password_data)
    background_tasks.add_task(send_password_has_been_reset_email, email_to=user.email, user_id=user.id)
    logger.info(f"User {user.email} has reset their password.")

    # log the user out by deleting any existing auth cookies.
    # TODO - when token blacklisting in place, blacklist the refresh tokens here.
    response = templates.TemplateResponse(
        request=request,
        name="auth_pages/login.html",
        context={"success": "Your password has been reset successfully, you can now log in."},
    )
    return delete_auth_cookies(response=response)
