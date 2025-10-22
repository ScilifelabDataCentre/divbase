"""
Global exception handlers for FastAPI.

Handles both API and frontend requests, returning JSON or HTML responses as appropriate.

The idea of centralising this is to:
1. Handle logging of exceptions all in one place.
2. Control what the users sees (not too much info and no accidental leakage)
3. Less work/duplication in the routes themselves, just raise the exception and the handler makes it pretty.

# TODO - see if the logging output is good.
TODO - generic 500 handler?
"""

import logging

from fastapi import FastAPI, Request, status
from fastapi.responses import JSONResponse, RedirectResponse

from divbase_api.exceptions import (
    AuthenticationError,
    AuthorizationError,
    BucketVersionAlreadyExistsError,
    BucketVersioningFileAlreadyExistsError,
    BucketVersionNotFoundError,
    ProjectCreationError,
    ProjectMemberNotFoundError,
    ProjectNotFoundError,
    TooManyObjectsInRequestError,
    UserRegistrationError,
)
from divbase_api.frontend_routes.core import templates

logger = logging.getLogger(__name__)


def is_api_request(request: Request) -> bool:
    """Helper function to check if request comes from frontend or API."""
    return request.url.path.startswith("/api/")


def render_error_page(request: Request, message: str, status_code: int = 500):
    """Helper function to render the generic error page for frontend requests."""
    return templates.TemplateResponse(
        request=request,
        name="error.html",
        context={
            "request": request,
            "error_message": message,
        },
        status_code=status_code,
    )


async def authentication_error_handler(request: Request, exc: AuthenticationError):
    logger.info(f"Authentication failed for {request.method} {request.url.path}: {exc.message}", exc_info=True)

    if is_api_request(request):
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": exc.message, "type": "authentication_error"},
            headers=exc.headers,
        )
    else:
        return RedirectResponse(url="/auth/login", status_code=status.HTTP_302_FOUND)


async def authorization_error_handler(request: Request, exc: AuthorizationError):
    logger.info(f"Authorization failed for {request.method} {request.url.path}: {exc.message}", exc_info=True)
    if is_api_request(request):
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": exc.message, "type": "authorization_error"},
            headers=exc.headers,
        )
    else:
        return RedirectResponse(url="/auth/login", status_code=status.HTTP_302_FOUND)


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
    logger.warning(f"Project not found for {request.method} {request.url.path}: {exc.message}", exc_info=True)
    if is_api_request(request):
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": exc.message, "type": "project_not_found_error"},
            headers=exc.headers,
        )
    else:
        return RedirectResponse(url="/projects", status_code=status.HTTP_302_FOUND)


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
        return render_error_page(request, exc.message, status_code=exc.status_code)


async def bucket_versioning_file_exists_error_handler(request: Request, exc: BucketVersioningFileAlreadyExistsError):
    logger.warning(
        f"Bucket versioning file already exists for {request.method} {request.url.path}: {exc.message}", exc_info=True
    )

    if is_api_request(request):
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": exc.message, "type": "bucket_versioning_file_already_exists_error"},
            headers=exc.headers,
        )
    else:
        return render_error_page(request, exc.message, status_code=exc.status_code)


async def bucket_version_exists_error_handler(request: Request, exc: BucketVersionAlreadyExistsError):
    logger.warning(
        f"Bucket version already exists for {request.method} {request.url.path}: {exc.message}", exc_info=True
    )

    if is_api_request(request):
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": exc.message, "type": "bucket_version_already_exists_error"},
            headers=exc.headers,
        )
    else:
        return render_error_page(request, exc.message, status_code=exc.status_code)


async def bucket_version_not_found_error_handler(request: Request, exc: BucketVersionNotFoundError):
    logger.warning(f"Bucket version not found for {request.method} {request.url.path}: {exc.message}", exc_info=True)

    if is_api_request(request):
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": exc.message, "type": "bucket_version_not_found_error"},
            headers=exc.headers,
        )
    else:
        return render_error_page(request, exc.message, status_code=exc.status_code)


def register_exception_handlers(app: FastAPI) -> None:
    """
    Register all exception handlers with FastAPI app.

    Type errors ignored (https://github.com/fastapi/fastapi/discussions/11741)
    """
    app.add_exception_handler(AuthenticationError, authentication_error_handler)  # type: ignore
    app.add_exception_handler(AuthorizationError, authorization_error_handler)  # type: ignore
    app.add_exception_handler(UserRegistrationError, user_registration_error_handler)  # type: ignore
    app.add_exception_handler(ProjectNotFoundError, project_not_found_error_handler)  # type: ignore
    app.add_exception_handler(ProjectMemberNotFoundError, project_member_not_found_error_handler)  # type: ignore
    app.add_exception_handler(ProjectCreationError, project_creation_error_handler)  # type: ignore
    app.add_exception_handler(TooManyObjectsInRequestError, too_many_objects_in_request_error_handler)  # type: ignore
    app.add_exception_handler(
        BucketVersioningFileAlreadyExistsError,
        bucket_versioning_file_exists_error_handler,  # type: ignore
    )
    app.add_exception_handler(BucketVersionAlreadyExistsError, bucket_version_exists_error_handler)  # type: ignore
    app.add_exception_handler(BucketVersionNotFoundError, bucket_version_not_found_error_handler)  # type: ignore
