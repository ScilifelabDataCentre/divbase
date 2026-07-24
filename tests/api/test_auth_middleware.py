"""Tests for middleware.py"""

from divbase_lib.divbase_constants import CLI_VERSION_HEADER_KEY


async def test_expected_headers_added_to_every_response(divbase_client):
    response = await divbase_client.get("/api/v1/core/health")
    assert response.status_code == 200
    assert response.headers["X-Request-ID"]
    assert response.headers["X-Content-Type-Options"] == "nosniff"


async def test_outdated_cli_version_is_rejected(divbase_client):
    """Rejected by CLIVersionMiddleware, before reaching the route handler."""
    response = await divbase_client.get("/api/v1/core/health", headers={CLI_VERSION_HEADER_KEY: "0.0.1"})
    assert response.status_code == 400
    assert response.json()["type"] == "CLIVersionOutdatedError"


async def test_request_with_disallowed_host_header_is_rejected(divbase_client):
    """Rejected by TrustedHostMiddleware, before reaching the route handler."""
    response = await divbase_client.get("/api/v1/core/health", headers={"Host": "evil.example.com"})
    assert response.status_code == 400
