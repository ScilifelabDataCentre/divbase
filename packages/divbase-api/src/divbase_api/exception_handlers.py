"""
Global exception handlers for FastAPI.

Handles both API and frontend requests, returning JSON or HTML responses as appropriate.

The idea of centralising this is to:
1. Handle logging of exceptions all in one place.
2. Control what the users sees (not too much info and no accidental leakage)
3. Less work/duplication in the routes themselves, just raise the exception and the handler makes it pretty.
"""

import logging

from fastapi import FastAPI, Request, status
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse, RedirectResponse
from starlette.exceptions import HTTPException

from divbase_api.db import get_db
from divbase_api.deps import _authenticate_frontend_user_from_tokens
from divbase_api.exceptions import (
    AuthenticationError,
    AuthorizationError,
    DownloadedFileChecksumMismatchError,
    ObjectDoesNotExistError,
    ProjectCreationError,
    ProjectMemberNotFoundError,
    ProjectNotFoundError,
    ProjectVersionAlreadyExistsError,
    ProjectVersionCreationError,
    ProjectVersionNotFoundError,
    TaskNotFoundInBackendError,
    TooManyObjectsInRequestError,
    UserRegistrationError,
    VCFDimensionsEntryMissingError,
)
from divbase_api.frontend_routes.core import templates
from divbase_api.models.users import UserDB

logger = logging.getLogger(__name__)


def is_api_request(request: Request) -> bool:
    """Helper function to check if request comes from frontend or API."""
    return request.url.path.startswith("/api/")


async def get_current_user_from_request_object(request: Request) -> UserDB | None:
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


async def render_error_page(
    request: Request,
    message: str,
    status_code: int = status.HTTP_500_INTERNAL_SERVER_ERROR,
):
    """Helper function to render the generic error page for frontend requests."""
    current_user = await get_current_user_from_request_object(request)
    return templates.TemplateResponse(
        request=request,
        name="error.html",
        context={
            "request": request,
            "error_message": message,
            "current_user": current_user,
        },
        status_code=status_code,
    )


async def global_exception_handler(request: Request, exc: Exception):
    """
    Handle unexpected exceptions globally. - in the ideal world this is never be triggered
    """
    logger.error(f"Unexpected Error occurred for: {request.method} {request.url.path}: {exc}", exc_info=True)
    return JSONResponse(
        status_code=500,
        content={
            "detail": "An unexpected error occurred. Please try again later.",
            "type": "server_error",
        },
    )


async def authentication_error_handler(request: Request, exc: AuthenticationError):
    logger.info(f"Authentication failed for {request.method} {request.url.path}: {exc.message}", exc_info=False)

    if is_api_request(request):
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": exc.message, "type": "authentication_error"},
            headers=exc.headers,
        )
    else:
        return RedirectResponse(url="/login", status_code=status.HTTP_302_FOUND)


async def authorization_error_handler(request: Request, exc: AuthorizationError):
    logger.info(f"Authorization failed for {request.method} {request.url.path}: {exc.message}", exc_info=False)
    if is_api_request(request):
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": exc.message, "type": "authorization_error"},
            headers=exc.headers,
        )
    else:
        return RedirectResponse(url="/login", status_code=status.HTTP_302_FOUND)


async def user_registration_error_handler(request: Request, exc: UserRegistrationError):
    logger.warning(f"User registration failed for {request.method} {request.url.path}: {exc.message}", exc_info=True)

    if is_api_request(request):
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": exc.user_message, "type": "user_registration_error"},
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


async def project_not_found_error_handler(request: Request, exc: ProjectNotFoundError):
    logger.info(f"Project not found for {request.method} {request.url.path}: {exc.message}", exc_info=False)
    if is_api_request(request):
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": exc.message, "type": "project_not_found_error"},
            headers=exc.headers,
        )
    else:
        return await render_error_page(
            request, "Project not found or you don't have access.", status_code=exc.status_code
        )


async def project_member_not_found_error_handler(request: Request, exc: ProjectMemberNotFoundError):
    logger.warning(f"Project member not found for {request.method} {request.url.path}: {exc.message}", exc_info=True)
    if is_api_request(request):
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": exc.message, "type": "project_member_not_found_error"},
            headers=exc.headers,
        )
    else:
        return RedirectResponse(url="/projects", status_code=status.HTTP_302_FOUND)


async def project_creation_error_handler(request: Request, exc: ProjectCreationError):
    logger.warning(f"Project creation failed for {request.method} {request.url.path}: {exc.message}", exc_info=True)
    if is_api_request(request):
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": exc.message, "type": "project_creation_error"},
            headers=exc.headers,
        )
    else:
        return RedirectResponse(url="/projects", status_code=status.HTTP_302_FOUND)


async def too_many_objects_in_request_error_handler(request: Request, exc: TooManyObjectsInRequestError):
    logger.warning(f"Too many objects in request for {request.method} {request.url.path}: {exc.message}", exc_info=True)

    if is_api_request(request):
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": exc.message, "type": "too_many_objects_in_request_error"},
            headers=exc.headers,
        )
    else:
        return await render_error_page(request, exc.message, status_code=exc.status_code)


async def project_version_creation_error_handler(request: Request, exc: ProjectVersionCreationError):
    logger.warning(
        f"Project version creation error for {request.method} {request.url.path}: {exc.message}", exc_info=True
    )

    if is_api_request(request):
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": exc.message, "type": "project_version_creation_error"},
            headers=exc.headers,
        )
    else:
        return await render_error_page(request, exc.message, status_code=exc.status_code)


async def project_version_already_exists_error_handler(request: Request, exc: ProjectVersionAlreadyExistsError):
    logger.info(
        f"Project version already exists error for {request.method} {request.url.path}: {exc.message}", exc_info=True
    )

    if is_api_request(request):
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": exc.message, "type": "project_version_already_exists_error"},
            headers=exc.headers,
        )
    else:
        return await render_error_page(request, exc.message, status_code=exc.status_code)


async def project_version_not_found_error_handler(request: Request, exc: ProjectVersionNotFoundError):
    logger.warning(
        f"Project version not found error for {request.method} {request.url.path}: {exc.message}", exc_info=True
    )

    if is_api_request(request):
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": exc.message, "type": "project_version_not_found_error"},
            headers=exc.headers,
        )
    else:
        return await render_error_page(request, exc.message, status_code=exc.status_code)


async def generic_http_exception_handler(request: Request, exc: HTTPException):
    """
    Generic handler for HTTP exceptions not caught by the specific handlers above.

    401 and 403s are handled by custom exceptions above.

    Note that we have to import starlette.exceptions.HTTPException here not FastAPI's HTTPException.
    """
    if is_api_request(request):
        if exc.status_code != 404:
            logger.error(
                f"HTTP {exc.status_code} error for {request.method} {request.url.path}: {exc.detail}", exc_info=True
            )
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": exc.detail, "type": "http_error"},
        )

    # (Frontend request)
    if exc.status_code == 404:
        current_user = await get_current_user_from_request_object(request)
        return templates.TemplateResponse(
            request=request,
            name="404.html",
            context={"request": request, "current_user": current_user},
            status_code=status.HTTP_404_NOT_FOUND,
        )
    else:
        logger.error(
            f"HTTP {exc.status_code} error for {request.method} {request.url.path}: {exc.detail}", exc_info=True
        )
        return await render_error_page(
            request, "An unexpected error occurred. Please try again later.", exc.status_code
        )


async def request_validation_error_handler(request: Request, exc: RequestValidationError):
    """When a request contains invalid data, FastAPI internally raises a RequestValidationError"""
    logger.error(f"Request validation error for {request.method} {request.url.path}: {exc.errors()}", exc_info=True)

    if is_api_request(request):
        return JSONResponse(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            content={"detail": exc.errors(), "type": "request_validation_error"},
        )
    else:
        return await render_error_page(
            request=request,
            message="Badly formatted request. Please check your input and try again.",
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
        )


async def vcf_dimensions_entry_missing_error_handler(request: Request, exc: VCFDimensionsEntryMissingError):
    logger.info(f"VCF dimensions entry missing for {request.method} {request.url.path}: {exc.message}", exc_info=False)

    if is_api_request(request):
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": exc.message, "type": "vcf_dimensions_entry_missing_error"},
            headers=exc.headers,
        )
    else:
        return await render_error_page(request, exc.message, status_code=exc.status_code)


async def task_not_found_in_backend_error_handler(request: Request, exc: TaskNotFoundInBackendError):
    logger.warning(
        f"Task ID not found in results backend for {request.method} {request.url.path}: {exc.message}", exc_info=True
    )
    if is_api_request(request):
        return JSONResponse(
            status_code=status.HTTP_410_GONE,
            content={"detail": exc.message, "type": "task_not_found_in_backend_error"},
            headers=exc.headers,
        )
    else:
        return await render_error_page(request, exc.message, status_code=status.HTTP_410_GONE)


async def downloaded_file_checksum_mismatch_error_handler(request: Request, exc: DownloadedFileChecksumMismatchError):
    logger.warning(
        f"Downloaded file checksum mismatch error for {request.method} {request.url.path}: {exc.message}", exc_info=True
    )
    if is_api_request(request):
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": exc.message, "type": "downloaded_file_checksum_mismatch_error"},
            headers=exc.headers,
        )
    else:
        return await render_error_page(request, exc.message, status_code=exc.status_code)


async def object_does_not_exist_error_handler(request: Request, exc: ObjectDoesNotExistError):
    logger.debug(f"Object does not exist error for {request.method} {request.url.path}: {exc.message}", exc_info=False)
    if is_api_request(request):
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": exc.message, "type": "object_does_not_exist_error"},
            headers=exc.headers,
        )
    else:
        return await render_error_page(
            request, "The requested file does not exist or you don't have access.", status_code=exc.status_code
        )


def register_exception_handlers(app: FastAPI) -> None:
    """
    Register all exception handlers with FastAPI app.

    Type errors ignored (https://github.com/fastapi/fastapi/discussions/11741)

    NOTE: error handlers need to be defined above this function, otherwise they will not work.
    """
    app.add_exception_handler(AuthenticationError, authentication_error_handler)  # type: ignore
    app.add_exception_handler(AuthorizationError, authorization_error_handler)  # type: ignore
    app.add_exception_handler(UserRegistrationError, user_registration_error_handler)  # type: ignore
    app.add_exception_handler(ProjectNotFoundError, project_not_found_error_handler)  # type: ignore
    app.add_exception_handler(ProjectMemberNotFoundError, project_member_not_found_error_handler)  # type: ignore
    app.add_exception_handler(ProjectCreationError, project_creation_error_handler)  # type: ignore
    app.add_exception_handler(ProjectVersionAlreadyExistsError, project_version_already_exists_error_handler)  # type: ignore
    app.add_exception_handler(TooManyObjectsInRequestError, too_many_objects_in_request_error_handler)  # type: ignore
    app.add_exception_handler(ProjectVersionCreationError, project_version_creation_error_handler)  # type: ignore
    app.add_exception_handler(ProjectVersionNotFoundError, project_version_not_found_error_handler)  # type: ignore
    app.add_exception_handler(RequestValidationError, request_validation_error_handler)  # type: ignore
    app.add_exception_handler(VCFDimensionsEntryMissingError, vcf_dimensions_entry_missing_error_handler)  # type: ignore
    app.add_exception_handler(TaskNotFoundInBackendError, task_not_found_in_backend_error_handler)  # type: ignore
    app.add_exception_handler(DownloadedFileChecksumMismatchError, downloaded_file_checksum_mismatch_error_handler)  # type: ignore
    app.add_exception_handler(ObjectDoesNotExistError, object_does_not_exist_error_handler)  # type: ignore

    # These cover more generic/unexpected HTTP errors - the exceptions above take precedence
    app.add_exception_handler(HTTPException, generic_http_exception_handler)  # type: ignore
    app.add_exception_handler(Exception, global_exception_handler)  # type: ignore
