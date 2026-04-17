"""
Playwright tests for the Personal Access Tokens (PATs) frontend pages.

Covers the PAT list page, token creation form, token revocation, and server-side validation errors.
"""

import re

import pytest
from playwright.sync_api import Browser, Page, expect
from sqlalchemy import delete, select

from divbase_api.crud.personal_access_tokens import PAT_MAX_ACTIVE_TOKENS
from divbase_api.models.personal_access_tokens import PersonalAccessTokenDB
from divbase_api.models.users import UserDB
from divbase_api.security import generate_personal_access_token, hash_personal_access_token

from .conftest import FRONTEND_BASE_URL, _create_logged_in_user_page, navigate_to

EDIT_USER_EMAIL = "edit@divbase.se"


@pytest.fixture
def logged_in_edit_user_pat_page(browser: Browser, EXISTING_ACCOUNTS):
    """Logged-in edit user page, isolated browser context for PAT tests."""
    page = _create_logged_in_user_page(browser, EXISTING_ACCOUNTS["EDIT_USER"])
    yield page
    page.close()


@pytest.fixture(autouse=True)
def clean_edit_user_pats(db_session_sync):
    """
    Delete all PATs for the edit user before and after each test to ensure test isolation.
    This prevents PAT names from colliding across test runs within the same session.
    """
    user = db_session_sync.execute(select(UserDB).where(UserDB.email == EDIT_USER_EMAIL)).scalar_one()
    db_session_sync.execute(delete(PersonalAccessTokenDB).where(PersonalAccessTokenDB.user_id == user.id))
    db_session_sync.commit()
    yield
    db_session_sync.execute(delete(PersonalAccessTokenDB).where(PersonalAccessTokenDB.user_id == user.id))
    db_session_sync.commit()


def fill_and_submit_new_pat_form(
    page: Page,
    name: str,
    expires_in_days: str = "30",
    task_history: bool = False,
) -> None:
    """
    Fill in and submit the new PAT creation form.

    Explicitly sets both scope checkboxes so the test result is
    deterministic regardless of the server-side default rendering.
    """
    navigate_to(page, "/pats/new")
    page.locator("#name").fill(name)
    page.locator("#expires_in_days").select_option(expires_in_days)

    if task_history:
        page.locator("#scope_task_history").check()
    else:
        page.locator("#scope_task_history").uncheck()

    page.get_by_role("button", name="Generate token").click()


def test_pat_list_page_renders_when_logged_in(logged_in_edit_user_pat_page: Page):
    """The PAT list page heading and 'Generate new token' link are visible when logged in."""
    navigate_to(logged_in_edit_user_pat_page, "/pats")
    expect(logged_in_edit_user_pat_page.get_by_role("heading", name="Personal Access Tokens")).to_be_visible()
    expect(logged_in_edit_user_pat_page.get_by_role("link", name=re.compile("Generate new token"))).to_be_visible()


def test_new_pat_form_renders(logged_in_edit_user_pat_page: Page):
    """The new PAT form shows the name field, expiration select, and submit button."""
    navigate_to(logged_in_edit_user_pat_page, "/pats/new")
    expect(logged_in_edit_user_pat_page.get_by_role("heading", name="New Personal Access Token")).to_be_visible()
    expect(logged_in_edit_user_pat_page.locator("#name")).to_be_visible()
    expect(logged_in_edit_user_pat_page.locator("#expires_in_days")).to_be_visible()
    expect(logged_in_edit_user_pat_page.get_by_role("button", name="Generate token")).to_be_visible()


def test_create_full_access_pat_shows_token(logged_in_edit_user_pat_page: Page):
    """
    Creating a PAT with no scope restrictions shows the 'Token generated' page,
    the copy-now warning, and the raw token (which must start with 'divbase_pat_').
    """
    fill_and_submit_new_pat_form(logged_in_edit_user_pat_page, name="my-hpc-token")

    expect(logged_in_edit_user_pat_page.get_by_role("heading", name="Token generated")).to_be_visible()
    expect(logged_in_edit_user_pat_page.get_by_text("Copy your token now, it will not be shown again.")).to_be_visible()
    # The raw token element contains the token value; verify its prefix
    token_text = logged_in_edit_user_pat_page.get_by_text(re.compile(r"divbase_pat_")).first.text_content()
    assert token_text is not None
    assert token_text.strip().startswith("divbase_pat_")


def test_created_pat_appears_in_list(logged_in_edit_user_pat_page: Page):
    """After creating a PAT and clicking 'Done', the token name appears in the list page."""
    fill_and_submit_new_pat_form(logged_in_edit_user_pat_page, name="list-check-token")

    logged_in_edit_user_pat_page.get_by_role("link", name=re.compile("Done")).click()
    expect(logged_in_edit_user_pat_page).to_have_url(f"{FRONTEND_BASE_URL}/pats/")
    # 1 for token name, 1 for hidden modal, if user clicks to delete the token
    expect(logged_in_edit_user_pat_page.get_by_text("list-check-token")).to_have_count(2)


def test_create_pat_with_never_expiry(logged_in_edit_user_pat_page: Page):
    """A PAT created with 'never' expiry shows 'Never' in the list."""
    fill_and_submit_new_pat_form(logged_in_edit_user_pat_page, name="no-expiry-token", expires_in_days="never")
    logged_in_edit_user_pat_page.get_by_role("link", name=re.compile("Done")).click()

    row = logged_in_edit_user_pat_page.get_by_role("row").filter(has_text="no-expiry-token")
    expect(row.get_by_text("Never").first).to_be_visible()


def test_create_pat_with_no_task_history_permissions(logged_in_edit_user_pat_page: Page):
    """A PAT created with no task_history permissions does not show the badge in the list."""
    fill_and_submit_new_pat_form(
        logged_in_edit_user_pat_page,
        name="scoped-token",
        task_history=False,
    )
    logged_in_edit_user_pat_page.get_by_role("link", name=re.compile("Done")).click()

    row = logged_in_edit_user_pat_page.get_by_role("row").filter(has_text="scoped-token")
    expect(row.get_by_text("task_history")).not_to_be_visible()


def test_revoke_pat_shows_success_message(logged_in_edit_user_pat_page: Page):
    """Revoking a PAT via the confirmation modal shows a success flash and removes the token from the list."""
    fill_and_submit_new_pat_form(logged_in_edit_user_pat_page, name="revoke-me-token")
    logged_in_edit_user_pat_page.get_by_role("link", name=re.compile("Done")).click()

    logged_in_edit_user_pat_page.get_by_role("button", name="Revoke").first.click()
    logged_in_edit_user_pat_page.get_by_role("button", name="Revoke token").click()

    expect(logged_in_edit_user_pat_page).to_have_url(f"{FRONTEND_BASE_URL}/pats/?success=Token+successfully+revoked")
    expect(logged_in_edit_user_pat_page.get_by_text(re.compile("successfully revoked", re.IGNORECASE))).to_be_visible()
    expect(logged_in_edit_user_pat_page.get_by_text("revoke-me-token")).not_to_be_visible()


def test_create_pat_duplicate_name_shows_error(logged_in_edit_user_pat_page: Page):
    """Submitting a PAT with a name already in use shows a duplicate name error."""
    fill_and_submit_new_pat_form(logged_in_edit_user_pat_page, name="duplicate-name")
    logged_in_edit_user_pat_page.get_by_role("link", name=re.compile("Done")).click()

    fill_and_submit_new_pat_form(logged_in_edit_user_pat_page, name="duplicate-name")
    expect(logged_in_edit_user_pat_page).to_have_url(f"{FRONTEND_BASE_URL}/pats/new")
    expect(
        logged_in_edit_user_pat_page.get_by_text(re.compile("already have an active token", re.IGNORECASE))
    ).to_be_visible()


def test_create_pat_limit_exceeded_shows_error(logged_in_edit_user_pat_page: Page, db_session_sync):
    user = db_session_sync.execute(select(UserDB).where(UserDB.email == EDIT_USER_EMAIL)).scalar_one()
    for i in range(PAT_MAX_ACTIVE_TOKENS):
        pat = PersonalAccessTokenDB(
            user_id=user.id,
            name=f"pre-filled-pat-{i}",
            hashed_token=hash_personal_access_token(generate_personal_access_token()),
            permissions={"all_projects": True, "projects": {}, "task_history": True},
        )
        db_session_sync.add(pat)
    db_session_sync.commit()

    fill_and_submit_new_pat_form(logged_in_edit_user_pat_page, name="one-too-many")
    expect(logged_in_edit_user_pat_page).to_have_url(f"{FRONTEND_BASE_URL}/pats/new")
    expect(
        logged_in_edit_user_pat_page.get_by_role("alert").filter(
            has_text=re.compile("reached the maximum", re.IGNORECASE)
        )
    ).to_be_visible()


def test_create_pat_with_role_exceeds_user_role_shows_error(logged_in_edit_user_pat_page: Page):
    """Selecting a project-specific role above the user's own role gives a validation error."""
    navigate_to(logged_in_edit_user_pat_page, "/pats/new")
    logged_in_edit_user_pat_page.locator("#name").fill("role-exceeds-test")
    logged_in_edit_user_pat_page.locator("#expires_in_days").select_option("30")

    # Switch to "specific projects" mode
    logged_in_edit_user_pat_page.locator('input[name="project_access_mode"][value="specific"]').click()

    # Select "manage" for the first project in the table (user only has edit)
    first_project_select = logged_in_edit_user_pat_page.locator("select[name^='project_']").first
    first_project_select.wait_for()
    first_project_select.select_option("manage")

    logged_in_edit_user_pat_page.get_by_role("button", name="Generate token").click()

    expect(logged_in_edit_user_pat_page).to_have_url(f"{FRONTEND_BASE_URL}/pats/new")
    expect(
        logged_in_edit_user_pat_page.get_by_text(re.compile("cannot exceed your membership role", re.IGNORECASE))
    ).to_be_visible()
