"""
Global exception handlers for FastAPI.

Handles both API and frontend requests, returning JSON or HTML responses as appropriate.

The idea of centralising this is to:
1. Handle logging of exceptions all in one place.
2. Control what the users sees (not too much info and no accidental leakage)
3. Less work/duplication in the routes themselves, just raise the exception and the handler makes it pretty.
"""

import logging

import structlog
from fastapi import FastAPI, Request, status
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse, RedirectResponse, Response
from starlette.exceptions import HTTPException

from divbase_api.db import get_db
from divbase_api.deps import _authenticate_frontend_user_from_tokens
from divbase_api.exceptions import DivBaseAPIException, UserRegistrationError
from divbase_api.frontend_routes.core import templates
from divbase_api.models.users import UserDB

logger = structlog.get_logger(__name__)


def _is_api_request(request: Request) -> bool:
    """Helper function to check if request comes from frontend or API."""
    return request.url.path.startswith("/api/")


def _log_exception(
    request: Request, exc: Exception, log_level: int, prefix: str, include_exc_info: bool = True
) -> None:
    """Helper function to log exception with request details"""
    event = f"{prefix} Details: {request.method} {request.url.path}: {exc}"
    logger.log(
        log_level,
        event,
        exc_info=include_exc_info,
    )


async def _get_current_user_from_request_object(request: Request) -> UserDB | None:
    """Helper function to get current user from request object"""
    access_token = request.cookies.get("access_token")
    refresh_token = request.cookies.get("refresh_token")
    if not access_token and not refresh_token:
        return None

    async for db in get_db():
        user = await _authenticate_frontend_user_from_tokens(
            access_token=access_token, refresh_token=refresh_token, db=db
        )
    return user


def _add_request_id_header(request: Request, response: Response) -> Response:
    """
    helper fn to add X-Request-ID to an existing response
    Needed in the global/generic exception handler as the middleware approach of appending the header will not work as won't be reached.
    """
    request_id = getattr(request.state, "request_id", None)
    if request_id:
        response.headers["X-Request-ID"] = request_id
    return response


async def _render_error_page(
    request: Request,
    message: str,
    error_template: str = "error.html",
    status_code: int = status.HTTP_500_INTERNAL_SERVER_ERROR,
):
    """Helper function to render an error page for frontend requests."""
    current_user = await _get_current_user_from_request_object(request)
    response = templates.TemplateResponse(
        request=request,
        name=error_template,
        context={
            "request": request,
            "error_message": message,
            "current_user": current_user,
        },
        status_code=status_code,
    )
    return _add_request_id_header(request, response)


async def divbase_api_exception_handler(request: Request, exc: DivBaseAPIException) -> Response:
    """
    The default exception handler for any DivBaseAPIException subclass.
    Each DivBaseAPIException child class declares how it should be handled using class attributes.

    Certain exceptions are handled individually as they have special requirements (e.g. UserRegistrationError).
    """
    _log_exception(
        request=request, exc=exc, log_level=exc.log_level, prefix=exc.error_type, include_exc_info=exc.include_exc_info
    )

    if _is_api_request(request):
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": exc.message, "type": exc.error_type},
            headers=exc.headers,
        )

    if exc.frontend_redirect_url:
        return RedirectResponse(url=exc.frontend_redirect_url, status_code=status.HTTP_302_FOUND)

    return await _render_error_page(
        request=request,
        message=exc.frontend_message or exc.message,
        status_code=exc.status_code,
        error_template="404.html" if exc.status_code == status.HTTP_404_NOT_FOUND else "error.html",
    )


async def unexpected_exception_handler(request: Request, exc: Exception):
    """
    Handle unexpected exceptions. In the ideal world this is never triggered.
    """
    prefix = f"Unexpected error '{type(exc).__name__}' occurred."
    _log_exception(request=request, exc=exc, log_level=logging.ERROR, prefix=prefix, include_exc_info=True)

    user_message = "An unexpected error occurred. Please try again later."
    if _is_api_request(request):
        return _add_request_id_header(
            request,
            JSONResponse(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                content={"detail": user_message, "type": "InternalServerError"},
            ),
        )
    else:
        return await _render_error_page(
            request=request, message=user_message, status_code=status.HTTP_500_INTERNAL_SERVER_ERROR
        )


async def user_registration_error_handler(request: Request, exc: UserRegistrationError):
    """
    UserRegistrationError has its own dedicated handler. As it needs to
    (1) show different user-facing message vs internal logging message and
    (2) can render 2 diff templates, depending on whether the error happened in the admin panel or public registration page.
    """
    prefix = f"{exc.error_type} occurred."
    _log_exception(
        request=request, exc=exc, log_level=exc.log_level, prefix=prefix, include_exc_info=exc.include_exc_info
    )

    if _is_api_request(request):
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": exc.user_message, "type": exc.error_type},
        )

    if request.url.path.startswith("/admin/"):
        # Admin user creation page
        return templates.TemplateResponse(
            request=request,
            name="admin_pages/create_user.html",
            context={
                "request": request,
                "error": exc.user_message,
            },
            status_code=exc.status_code,
        )

    # Public registration page
    return templates.TemplateResponse(
        request=request,
        name="auth_pages/register.html",
        context={
            "request": request,
            "error": exc.user_message,
        },
        status_code=exc.status_code,
    )


async def generic_http_exception_handler(request: Request, exc: HTTPException):
    """
    Generic handler for HTTP exceptions. 401 and 403s are handled by custom exceptions which inherit from DivBaseAPIException.

    NOTE: that we have to import starlette.exceptions.HTTPException here not FastAPI's HTTPException.
    """
    error_type = type(exc).__name__
    _log_exception(
        request=request,
        exc=exc,
        log_level=logging.ERROR if exc.status_code != 404 else logging.INFO,
        prefix=error_type,
        include_exc_info=exc.status_code != 404,
    )

    if _is_api_request(request):
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": exc.detail, "type": error_type},
            headers=exc.headers,
        )

    # (Frontend request)
    if exc.status_code == 404:
        current_user = await _get_current_user_from_request_object(request)
        return templates.TemplateResponse(
            request=request,
            name="404.html",
            context={"request": request, "current_user": current_user},
            status_code=status.HTTP_404_NOT_FOUND,
        )
    else:
        return await _render_error_page(
            request=request,
            message="An unexpected error occurred. Please try again later.",
            status_code=exc.status_code,
        )


async def request_validation_error_handler(request: Request, exc: RequestValidationError):
    """
    When a request contains invalid data, FastAPI internally raises a RequestValidationError

    NOTE: Pydantic returns a list of validation errors, so this handler works differently than others.
    """
    error_type = "RequestValidationError"
    logger.info(f"{error_type} for {request.method} {request.url.path}: {exc.errors()}", exc_info=False)

    # gives the user back multiple error messages if there is more than 1 validation error
    errors = [error["msg"] for error in exc.errors() if "msg" in error]
    if _is_api_request(request):
        return JSONResponse(
            status_code=status.HTTP_422_UNPROCESSABLE_CONTENT,
            content={"detail": errors, "type": error_type},
        )
    else:
        return await _render_error_page(
            request=request,
            message=f"Badly formatted request. Please check your input and try again. Details: {errors}",
            status_code=status.HTTP_422_UNPROCESSABLE_CONTENT,
        )


def register_exception_handlers(app: FastAPI) -> None:
    """
    Register all exception handlers with FastAPI app.

    Type errors ignored (https://github.com/fastapi/fastapi/discussions/11741)

    NOTE: error handlers need to be defined above this function, otherwise they will not work.
    NOTE: UserRegistrationError has its own dedicated handler (custom templates/messages) rather than
    going through divbase_api_exception_handler. (FastAPI chooses whichever registered exception handler is most specific.)
    """
    app.add_exception_handler(DivBaseAPIException, divbase_api_exception_handler)  # type: ignore
    app.add_exception_handler(UserRegistrationError, user_registration_error_handler)  # type: ignore
    app.add_exception_handler(RequestValidationError, request_validation_error_handler)  # type: ignore

    # These cover more generic/unexpected HTTP errors - the exceptions above take precedence
    app.add_exception_handler(HTTPException, generic_http_exception_handler)  # type: ignore
    app.add_exception_handler(Exception, unexpected_exception_handler)  # type: ignore
