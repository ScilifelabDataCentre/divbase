"""
Retry logic for the divbase-cli package.

Defines functions to determine whether to retry based on the type of exceptions raised.
Functions are used in stamina (package) decorators.
"""

import httpx

from divbase_cli.cli_exceptions import DivBaseAPIConnectionError, DivBaseAPIError


def retry_only_on_retryable_http_errors(exc: Exception) -> bool:
    """
    Used by stamina's (library for retries) decorators to determine whether to retry the function or not.
    We avoid retrying on HTTPStatusError for 4xx errors as no point (e.g. 404 Not Found or 403 Forbidden etc...).
    """
    if isinstance(exc, httpx.HTTPStatusError):
        return exc.response.status_code >= 500

    # Want to retry on other HTTPError (parent of HTTPStatusError),
    # as this includes timeouts, connection errors, etc.
    return isinstance(exc, httpx.HTTPError)


def retry_only_on_retryable_divbase_api_errors(exception: Exception) -> bool:
    """
    Retry condition function for stamina to only retry on retryable DivBase API errors.
    """
    if isinstance(exception, DivBaseAPIConnectionError):
        return True

    # Retry only for server errors (5xx) and rate limiting (429)
    return isinstance(exception, DivBaseAPIError) and (exception.status_code == 429 or exception.status_code >= 500)
