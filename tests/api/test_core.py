"""API tests for routes/core.py."""

from datetime import datetime, timedelta, timezone

from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.models.announcements import AnnouncementDB, AnnouncementLevel, AnnouncementTarget
from divbase_lib import __version__ as divbase_server_version
from divbase_lib.api_schemas.announcements import AnnouncementResponse
from divbase_lib.divbase_constants import CLI_VERSION_HEADER_KEY

ANNOUNCEMENTS_ENDPOINT = "/api/v1/core/announcements"


async def _create_announcement(
    db_session: AsyncSession,
    target: AnnouncementTarget,
    heading: str = "Test announcement",
    auto_expire_at: datetime | None = None,
) -> None:
    announcement = AnnouncementDB(
        heading=heading,
        message="Some details.",
        target=target,
        level=AnnouncementLevel.INFO,
        auto_expire_at=auto_expire_at,
    )
    db_session.add(announcement)
    await db_session.commit()


async def test_announcements_with_none_active_returns_empty_list(divbase_client):
    response = await divbase_client.get(ANNOUNCEMENTS_ENDPOINT)
    assert response.status_code == 200
    assert response.json() == []


async def test_announcements_returns_active_cli_and_both_targeted_announcements(divbase_client, db_session):
    """This endpoint is used by the CLI, web is for display on homepage"""
    await _create_announcement(db_session, heading="CLI only", target=AnnouncementTarget.CLI)
    await _create_announcement(db_session, heading="Both", target=AnnouncementTarget.BOTH)
    await _create_announcement(db_session, heading="Web only", target=AnnouncementTarget.WEB)

    response = await divbase_client.get(ANNOUNCEMENTS_ENDPOINT)
    assert response.status_code == 200
    announcements = [AnnouncementResponse.model_validate(a) for a in response.json()]
    headings = {a.heading for a in announcements}
    assert headings == {"CLI only", "Both"}


async def test_announcements_excludes_expired_announcement(divbase_client, db_session):
    """Test that expired announcements are not returned by the API."""
    await _create_announcement(
        db_session=db_session,
        target=AnnouncementTarget.CLI,
        heading="Expired",
        auto_expire_at=datetime.now(timezone.utc) - timedelta(days=1),
    )
    await _create_announcement(
        db_session=db_session,
        target=AnnouncementTarget.CLI,
        heading="Still active",
        auto_expire_at=datetime.now(timezone.utc) + timedelta(days=1),
    )

    response = await divbase_client.get(ANNOUNCEMENTS_ENDPOINT)
    assert response.status_code == 200
    announcements = [AnnouncementResponse.model_validate(a) for a in response.json()]
    headings = {a.heading for a in announcements}
    assert headings == {"Still active"}


async def test_announcements_adds_new_cli_version_notice_for_outdated_cli(divbase_client, monkeypatch):
    """Test that the API returns a info level notice that the user's CLI version can be upgraded."""
    # Ensure the middleware to reject outdated CLI versions does not interfere with this test (of the route handler).
    monkeypatch.setattr("divbase_api.middleware.cli_version_outdated", lambda cli_version: False)

    response = await divbase_client.get(ANNOUNCEMENTS_ENDPOINT, headers={CLI_VERSION_HEADER_KEY: "0.0.1"})
    assert response.status_code == 200
    announcements = [AnnouncementResponse.model_validate(a) for a in response.json()]
    assert len(announcements) == 1
    assert "A new version" in announcements[0].heading
    assert announcements[0].level == AnnouncementLevel.INFO


async def test_announcements_has_no_new_cli_version_notice_for_up_to_date_cli(divbase_client):
    """Test that the API does not return a notice for up-to-date or newer CLI versions."""
    # up-to-date CLI version
    response = await divbase_client.get(
        ANNOUNCEMENTS_ENDPOINT, headers={CLI_VERSION_HEADER_KEY: divbase_server_version}
    )
    assert response.status_code == 200
    assert response.json() == []
    # newer CLI version than server
    response = await divbase_client.get(ANNOUNCEMENTS_ENDPOINT, headers={CLI_VERSION_HEADER_KEY: "999.0.0"})
    assert response.status_code == 200
    assert response.json() == []
