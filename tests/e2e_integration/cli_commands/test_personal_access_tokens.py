"""
E2E tests validating that CLI commands work correctly when authenticated via Personal Access Tokens (PATs).

All tests use the `pat_factory` fixture to create a single PAT directly in the DB, which is auto-clean up after thetest.
and `monkeypatch` to set DIVBASE_API_PAT on cli_settings without touching the JWT-based test setup.

### Why only a subset of commands?
Almost all API endpoints are scoped based on project level permissions which is handled by the dependancy:
get_project_member in deps.py. Other more special cases like auth whoami and task-history user have their logic.

So below we cover a representative cmd for project scope (files ls) and then the special cases.
"""

from datetime import datetime, timedelta, timezone

import httpx
import pytest
from sqlalchemy import delete, select
from typer.testing import CliRunner

from divbase_api.models.personal_access_tokens import PersonalAccessTokenDB
from divbase_api.models.users import UserDB
from divbase_api.security import generate_personal_access_token, hash_personal_access_token
from divbase_cli.cli_config import cli_settings
from divbase_cli.cli_exceptions import DivBaseAPIError
from divbase_cli.divbase_cli import app

runner = CliRunner()

READ_USER_EMAIL = "read@divbase.se"
EDIT_USER_EMAIL = "edit@divbase.se"
MANAGE_USER_EMAIL = "manage@divbase.se"


FULL_ACCESS_PAT_PERMISSIONS = {
    "all_projects": True,
    "projects": {},
    "task_history": True,
}

NO_SCOPE_PERMISSIONS = {
    "all_projects": False,
    "projects": {},
    "task_history": False,
}
TASK_HISTORY_ONLY_PERMISSIONS = {
    "all_projects": False,
    "projects": {},
    "task_history": True,
}


def assert_401_error(result):
    """Helper to assert that a CLI result is an 401 error."""
    assert result.exit_code != 0
    assert isinstance(result.exception, DivBaseAPIError)
    assert "401" in str(result.exception)
    assert "authentication_error" in str(result.exception)


def assert_403_error(result):
    """Helper to assert that a CLI result is an 403 error."""
    assert result.exit_code != 0
    assert isinstance(result.exception, DivBaseAPIError)
    assert "403" in str(result.exception)
    assert "authorization_error" in str(result.exception)


@pytest.fixture
def pat_factory(db_session_sync):
    """Creates a single PAT directly in the DB and deletes it after the test."""
    pat_id: int | None = None

    def create_pat(user_email, permissions=None, expires_at=None):
        nonlocal pat_id
        user = db_session_sync.execute(select(UserDB).where(UserDB.email == user_email)).scalar_one()
        raw = generate_personal_access_token()
        pat = PersonalAccessTokenDB(
            user_id=user.id,
            name="test-pat",
            hashed_token=hash_personal_access_token(raw),
            permissions=permissions,
            expires_at=expires_at,
        )
        db_session_sync.add(pat)
        db_session_sync.commit()
        db_session_sync.refresh(pat)
        pat_id = pat.id
        return raw.get_secret_value()

    yield create_pat

    if pat_id is not None:
        db_session_sync.execute(delete(PersonalAccessTokenDB).where(PersonalAccessTokenDB.id == pat_id))
        db_session_sync.commit()


def test_full_access_pat_works_for_all_cmds(CONSTANTS, logged_out_user_with_existing_config, pat_factory, monkeypatch):
    """A full access PAT should work for all cms, see top for why this selection of cmds is enough."""
    raw_token = pat_factory(user_email=MANAGE_USER_EMAIL)
    monkeypatch.setattr(cli_settings, "DIVBASE_API_PAT", raw_token)

    result = runner.invoke(app, "auth whoami")
    assert result.exit_code == 0
    assert MANAGE_USER_EMAIL in result.output

    result = runner.invoke(app, "task-history user")
    assert result.exit_code == 0

    result = runner.invoke(app, f"task-history project {CONSTANTS['DEFAULT_PROJECT']}")
    assert result.exit_code == 0

    # don't have to specify the project as default will be taken from config
    result = runner.invoke(app, "task-history project")
    assert result.exit_code == 0

    result = runner.invoke(app, f"files ls --project {CONSTANTS['DEFAULT_PROJECT']}")
    assert result.exit_code == 0


def test_no_scope_pat_can_only_run_whoami(CONSTANTS, logged_out_user_with_existing_config, pat_factory, monkeypatch):
    """whoami is always accessible to any valid PAT regardless of scopes."""
    raw_token = pat_factory(
        user_email=EDIT_USER_EMAIL,
        permissions=NO_SCOPE_PERMISSIONS,
    )
    monkeypatch.setattr(cli_settings, "DIVBASE_API_PAT", raw_token)

    result = runner.invoke(app, "auth whoami")
    assert result.exit_code == 0
    assert EDIT_USER_EMAIL in result.output

    result = runner.invoke(app, "task-history user")
    assert_403_error(result)

    result = runner.invoke(app, f"files ls --project {CONSTANTS['DEFAULT_PROJECT']}")
    assert_403_error(result)


def test_task_history_only_can_only_run_task_history(
    CONSTANTS, logged_out_user_with_existing_config, pat_factory, monkeypatch
):
    raw_token = pat_factory(
        user_email=MANAGE_USER_EMAIL,
        permissions=TASK_HISTORY_ONLY_PERMISSIONS,
    )
    monkeypatch.setattr(cli_settings, "DIVBASE_API_PAT", raw_token)

    result = runner.invoke(app, "task-history user")
    assert result.exit_code == 0

    # (project-scoped task history cmd)
    result = runner.invoke(app, f"task-history project {CONSTANTS['DEFAULT_PROJECT']}")
    assert_403_error(result)

    result = runner.invoke(app, "auth whoami")
    assert result.exit_code == 0  # whoami is always accessible regardless of PAT scopes

    result = runner.invoke(app, f"files ls --project {CONSTANTS['DEFAULT_PROJECT']}")
    assert_403_error(result)


def test_project_scoped_pat(CONSTANTS, logged_out_user_with_existing_config, pat_factory, project_map, monkeypatch):
    """Validate works on project included in permissions and is rejected on excluded project."""
    included_project = CONSTANTS["DEFAULT_PROJECT"]
    excluded_project = CONSTANTS["NON_DEFAULT_PROJECT"]
    permissions = {
        "all_projects": False,
        "projects": {str(project_map[included_project]): "edit"},
        "task_history": False,
    }
    raw_token = pat_factory(user_email=EDIT_USER_EMAIL, permissions=permissions)
    monkeypatch.setattr(cli_settings, "DIVBASE_API_PAT", raw_token)

    result = runner.invoke(app, f"files ls --project {included_project}")
    assert result.exit_code == 0

    result = runner.invoke(app, f"files ls --project {excluded_project}")
    assert_403_error(result)


def test_project_scoped_task_history_included_project_succeeds(
    CONSTANTS, logged_out_user_with_existing_config, pat_factory, project_map, monkeypatch
):
    project = CONSTANTS["DEFAULT_PROJECT"]
    permissions = {
        "all_projects": False,
        "projects": {str(project_map[project]): "manage"},
        "task_history": True,
    }
    raw_token = pat_factory(user_email=MANAGE_USER_EMAIL, permissions=permissions)
    monkeypatch.setattr(cli_settings, "DIVBASE_API_PAT", raw_token)

    result = runner.invoke(app, f"task-history project {project}")
    assert result.exit_code == 0


def test_project_scoped_task_history_excluded_project_rejected(
    CONSTANTS, logged_out_user_with_existing_config, pat_factory, project_map, monkeypatch
):
    included_project = CONSTANTS["DEFAULT_PROJECT"]
    excluded_project = CONSTANTS["NON_DEFAULT_PROJECT"]
    permissions = {
        "all_projects": False,
        "projects": {str(project_map[included_project]): "edit"},
        "task_history": True,
    }
    raw_token = pat_factory(user_email=EDIT_USER_EMAIL, permissions=permissions)
    monkeypatch.setattr(cli_settings, "DIVBASE_API_PAT", raw_token)

    result = runner.invoke(app, f"task-history project {excluded_project}")
    assert_403_error(result)


def test_read_role_pat_rejects_edit_level_action(
    CONSTANTS, logged_out_user_with_existing_config, pat_factory, project_map, monkeypatch
):
    """A PAT capped at read role cannot submit a query (requires edit)."""
    project = CONSTANTS["DEFAULT_PROJECT"]
    permissions = {
        "all_projects": False,
        "projects": {str(project_map[project]): "read"},
        "task_history": False,
    }
    raw_token = pat_factory(user_email=EDIT_USER_EMAIL, permissions=permissions)
    monkeypatch.setattr(cli_settings, "DIVBASE_API_PAT", raw_token)

    result = runner.invoke(app, f"query tsv 'Area:West of Ireland' --project {project}")
    assert_403_error(result)


def test_edit_pat_cant_do_edit_actions_if_user_is_read_only(
    CONSTANTS, logged_out_user_with_existing_config, pat_factory, project_map, monkeypatch
):
    """
    Even if a PAT has edit permissions, if the user only has read permissions, edit actions should be rejected.

    NOTE: The frontend prevents creating a PAT with greater permissions than user's, but
    this is still possible if e.g. user's status is later downgraded after already creating a PAT.
    """
    project = CONSTANTS["DEFAULT_PROJECT"]
    permissions = {
        "all_projects": False,
        "projects": {str(project_map[project]): "edit"},
        "task_history": False,
    }
    raw_token = pat_factory(user_email=READ_USER_EMAIL, permissions=permissions)
    monkeypatch.setattr(cli_settings, "DIVBASE_API_PAT", raw_token)

    result = runner.invoke(app, f"query tsv 'Area:West of Ireland' --project {project}")
    assert_403_error(result)


def test_expired_pat_rejected(CONSTANTS, logged_out_user_with_existing_config, pat_factory, monkeypatch):
    expiry = datetime.now(timezone.utc) - timedelta(days=1)
    raw_token = pat_factory(user_email=EDIT_USER_EMAIL, expires_at=expiry)
    monkeypatch.setattr(cli_settings, "DIVBASE_API_PAT", raw_token)

    result = runner.invoke(app, "auth whoami")
    assert_401_error(result)

    result = runner.invoke(app, f"files ls --project {CONSTANTS['DEFAULT_PROJECT']}")
    assert_401_error(result)


def test_pat_does_not_work_on_frontend_endpoints(
    CONSTANTS, logged_out_user_with_existing_config, pat_factory, monkeypatch
):
    """PATs should not work on frontend endpoints"""
    raw_token = pat_factory(user_email=EDIT_USER_EMAIL, permissions=FULL_ACCESS_PAT_PERMISSIONS)

    home_url = "http://localhost:8001/"
    # technically don't need the token here
    response = httpx.get(home_url, headers={"Authorization": f"Bearer {raw_token}"})
    assert response.status_code == 200

    profile_url = "http://localhost:8001/profile/"
    response = httpx.get(profile_url, headers={"Authorization": f"Bearer {raw_token}"})
    assert response.status_code == 302  # redirected to login, so rejected as expected
