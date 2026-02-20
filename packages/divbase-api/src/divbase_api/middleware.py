"""
Middleware for Divbase API.

We do not need CORS middleware as the frontend has the same origin.

NOTE: The following are not implemented here as they are handled by the ingress in front of the API:
- Rate limiting
- ip address white listing for admin routes.
- Autoredirect from HTTP to HTTPS (so no need for HTTPSRedirectMiddleware)
"""

from urllib.parse import urlparse

from fastapi import FastAPI
from fastapi.middleware.gzip import GZipMiddleware
from starlette.middleware.base import BaseHTTPMiddleware, RequestResponseEndpoint
from starlette.middleware.trustedhost import TrustedHostMiddleware
from starlette.requests import Request
from starlette.responses import Response

from divbase_api.api_config import settings

IMAGE_FONT_EXTENSIONS = [".webp", ".svg", ".jpg", ".jpeg", ".png", ".woff", ".woff2", ".ttf"]
CSS_JS_EXTENSIONS = [".css", ".js"]


class CustomHeaderMiddleware(BaseHTTPMiddleware):
    async def dispatch(self, request: Request, call_next: RequestResponseEndpoint) -> Response:
        """
        Adds extra headers to all responses.
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


def register_middleware(app: FastAPI) -> None:
    """
    Register all middleware to the FastAPI app.

    NOTE: When modyfying/extending this function remember that Middleware order is important.

    So for the request-response cycle:
    request -> last added middleware -> ... -> first added middleware -> route handler
    route handler -> first added middleware -> ... -> last added middleware -> response
    """
    if settings.api.environment in ["local_dev", "test"]:
        allowed_hosts = ["localhost"]
    else:
        frontend_host = urlparse(settings.api.frontend_base_url).hostname
        if not frontend_host:
            raise ValueError("Could not parse hostname from FRONTEND_BASE_URL.")
        allowed_hosts = [frontend_host]

    app.add_middleware(middleware_class=GZipMiddleware, minimum_size=1000, compresslevel=5)
    app.add_middleware(CustomHeaderMiddleware)
    app.add_middleware(middleware_class=TrustedHostMiddleware, allowed_hosts=allowed_hosts)
