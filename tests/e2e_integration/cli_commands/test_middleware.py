"""
E2E tests for API middleware behavior.

These tests intentionally use real HTTP calls against the running docker test stack,
with no patching/mocking, to validate middleware behavior end-to-end.
"""

from concurrent.futures import ThreadPoolExecutor
from uuid import UUID

import httpx

from divbase_cli import __version__ as cli_version
from divbase_cli.cli_config import cli_settings
from divbase_lib.divbase_constants import CLI_VERSION_HEADER_KEY

HEALTH_URL = f"{cli_settings.DIVBASE_API_URL}/v1/core/health"


def _assert_valid_request_id(request_id: str) -> str:
    try:
        UUID(request_id)
    except ValueError as e:
        raise AssertionError(f"Invalid request ID, got: {request_id}. Expected a UUID string.") from e
    return request_id


def test_concurrent_health_requests_have_distinct_request_ids():
    """Three overlapping requests should each receive their own request id."""
    with ThreadPoolExecutor(max_workers=3) as executor:
        responses = list(executor.map(lambda _: httpx.get(HEALTH_URL, timeout=10), range(3)))

    request_ids: list[str] = []
    for response in responses:
        assert response.status_code == 200
        assert response.json() == {"status": "ok"}
        request_id = response.headers["X-Request-ID"]
        _assert_valid_request_id(request_id)
        request_ids.append(request_id)

    assert len(set(request_ids)) == 3


def test_cli_version_middleware_rejects_outdated_cli_header():
    """
    Validate outdated CLI version header is rejected.
    Should provide back a request ID header for both succesful and unsuccessful responses, and they should be different.
    """
    outdated_response = httpx.get(HEALTH_URL, headers={CLI_VERSION_HEADER_KEY: "0.0.0"}, timeout=10)
    assert outdated_response.status_code == 400
    assert outdated_response.json()["type"] == "cli_version_outdated_error"
    outdated_request_id = outdated_response.headers["X-Request-ID"]
    _assert_valid_request_id(outdated_request_id)

    ok_response = httpx.get(HEALTH_URL, headers={CLI_VERSION_HEADER_KEY: cli_version}, timeout=10)
    assert ok_response.status_code == 200
    assert ok_response.json() == {"status": "ok"}
    ok_request_id = ok_response.headers["X-Request-ID"]
    _assert_valid_request_id(ok_request_id)

    assert outdated_request_id != ok_request_id


def test_cli_version_header_is_optional():
    """
    Validate that:
     - Missing CLI version header should be accepted.
     - Malformed CLI version header should be rejected.
    """
    missing_header_response = httpx.get(HEALTH_URL, timeout=10)
    assert missing_header_response.status_code == 200
    assert missing_header_response.json() == {"status": "ok"}
    missing_request_id = missing_header_response.headers["X-Request-ID"]
    _assert_valid_request_id(missing_request_id)

    malformed_header_response = httpx.get(
        HEALTH_URL, headers={CLI_VERSION_HEADER_KEY: "not.a.valid.version"}, timeout=10
    )
    assert malformed_header_response.status_code == 400
    assert malformed_header_response.json()["type"] == "cli_version_outdated_error"
    malformed_request_id = malformed_header_response.headers["X-Request-ID"]
    _assert_valid_request_id(malformed_request_id)


def test_trusted_host_middleware_rejects_invalid_host():
    """Requests with an invalid Host header should be rejected by TrustedHostMiddleware."""
    response = httpx.get(HEALTH_URL, headers={"Host": "invalid.host.example"}, timeout=10)

    assert response.status_code == 400
    assert "Invalid host header" in response.text
