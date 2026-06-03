"""
Middleware for Divbase API.

We do not need CORS middleware as the frontend has the same origin.

NOTE: The following are not implemented here as they are handled by the ingress in front of the API:
- Rate limiting
- ip address white listing for admin routes.
- Autoredirect from HTTP to HTTPS (so no need for HTTPSRedirectMiddleware)
"""

import uuid
from urllib.parse import urlparse

import structlog
from fastapi import FastAPI
from fastapi.middleware.gzip import GZipMiddleware
from starlette.middleware.base import BaseHTTPMiddleware, RequestResponseEndpoint
from starlette.middleware.trustedhost import TrustedHostMiddleware
from starlette.requests import Request
from starlette.responses import JSONResponse, Response

from divbase_api.api_config import api_settings
from divbase_api.services.validate_cli_versions import cli_version_outdated
from divbase_lib.divbase_constants import CLI_VERSION_HEADER_KEY, LOCAL_DEV_ENVIRONMENTS

IMAGE_FONT_EXTENSIONS = [".webp", ".svg", ".jpg", ".jpeg", ".png", ".woff", ".woff2", ".ttf"]
CSS_JS_EXTENSIONS = [".css", ".js"]


class RequestContextMiddleware(BaseHTTPMiddleware):
    async def dispatch(self, request: Request, call_next: RequestResponseEndpoint) -> Response:
        """
        To every request that comes in, we bind(aka add) context into structlog's contextvars so every log line
        produced during the request-response cycle automatically includes these fields.

        2 benefits:
        - correlate every log message to the request.
        - user can give us the "X-Request-ID" in error report and we can search logs more easily.
        """
        request_id = str(uuid.uuid4())
        request.state.request_id = request_id
        structlog.contextvars.clear_contextvars()
        structlog.contextvars.bind_contextvars(request_id=request_id)
        response = await call_next(request)
        response.headers["X-Request-ID"] = request_id
        return response


class CustomHeaderMiddleware(BaseHTTPMiddleware):
    async def dispatch(self, request: Request, call_next: RequestResponseEndpoint) -> Response:
        """
        Adds extra headers to all responses after being processed by the route handler.
        Note that some other headers are added via the k8s ingress sitting in front of the API.
        """
        response = await call_next(request)
        path = request.url.path

        response.headers["X-Content-Type-Options"] = "nosniff"

        # everything below is only relevant for frontend routes
        if path.startswith("/api/"):
            return response

        response.headers["Content-Security-Policy"] = "frame-ancestors 'self'"
        response.headers["Referrer-Policy"] = "strict-origin-when-cross-origin"

        if any(path.endswith(ext) for ext in IMAGE_FONT_EXTENSIONS):
            response.headers["Cache-Control"] = "public, max-age=2592000"  # 30 days
        elif any(path.endswith(ext) for ext in CSS_JS_EXTENSIONS):
            # "no-cache" does not mean what you might think it means,
            # we are still caching, but browser checks with server if the file has changed first (via etag).
            # We don't do fingerprinting/versioning of static assets, so this makes sense.
            response.headers["Cache-Control"] = "no-cache"

        return response


class CLIVersionMiddleware(BaseHTTPMiddleware):
    async def dispatch(self, request: Request, call_next: RequestResponseEndpoint) -> Response:
        """
        Middleware to validate the CLI version from the request header.

        Will reject request if user cli version is too outdated before hitting the route handler.
        """
        cli_version = request.headers.get(CLI_VERSION_HEADER_KEY)

        if cli_version and cli_version_outdated(cli_version=cli_version):
            # in middleware, you can't rely on raising exceptions and passing to exception_handlers.py
            # You have to return a direct response instead.
            message = (
                "Your install of divbase is too outdated and no longer compatible with DivBase Server. "
                "You must first update your install of divbase in order to run any more commands. "
                "If you're not sure how to do that, you can find instructions on how to upgrade here: "
                f"{api_settings.general.mkdocs_site_url}/user-guides/installation"
            )
            body = {"detail": message, "type": "cli_version_outdated_error"}
            return JSONResponse(content=body, status_code=400)

        response = await call_next(request)
        return response


def register_middleware(app: FastAPI) -> None:
    """
    Register all middleware to the FastAPI app.

    NOTE: When modifying/extending this function remember that Middleware order is important.

    So for the request-response cycle:
    request -> last added middleware -> ... -> first added middleware -> route handler
    route handler -> first added middleware -> ... -> last added middleware -> response
    """
    if api_settings.general.environment in LOCAL_DEV_ENVIRONMENTS:
        allowed_hosts = ["localhost"]
    else:
        frontend_host = urlparse(api_settings.general.frontend_base_url).hostname
        if not frontend_host:
            raise ValueError("Could not parse hostname from FRONTEND_BASE_URL.")
        allowed_hosts = [frontend_host]

    app.add_middleware(middleware_class=GZipMiddleware, minimum_size=1000, compresslevel=5)
    app.add_middleware(CustomHeaderMiddleware)
    app.add_middleware(CLIVersionMiddleware)
    app.add_middleware(RequestContextMiddleware)
    app.add_middleware(middleware_class=TrustedHostMiddleware, allowed_hosts=allowed_hosts)
