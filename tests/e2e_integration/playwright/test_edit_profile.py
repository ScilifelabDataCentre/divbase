"""
Test edit profile page functionality using Playwright.
"""

import pytest
from playwright.sync_api import Page, expect

from divbase_api.api_constants import KNOWN_JOB_ROLES, SWEDISH_UNIVERSITIES

from .conftest import FRONTEND_BASE_URL, navigate_to


def _edit_profile(
    page: Page,
    new_name: str,
    new_organisation: str,
    new_role: str,
):
    """
    Helper function to edit the profile of a user.
    Handles the "Other" option for organisation and role.
    """
    navigate_to(page=page, path="/profile/edit")

    organisation_input = page.get_by_label("Organisation", exact=True)
    role_input = page.get_by_label("Role", exact=True)
    other_organisation_input = page.get_by_role("textbox", name="Your organisation")
    other_role_input = page.get_by_role("textbox", name="Your Role")

    page.get_by_role("textbox", name="Full Name").fill(new_name)

    if new_organisation in SWEDISH_UNIVERSITIES:
        organisation_input.select_option(new_organisation)
        expect(other_organisation_input).to_be_hidden()
    else:
        organisation_input.select_option("Other")
        expect(other_organisation_input).to_be_visible()
        other_organisation_input.fill(new_organisation)

    if new_role in KNOWN_JOB_ROLES:
        role_input.select_option(new_role)
        expect(other_role_input).to_be_hidden()
    else:
        role_input.select_option("Other")
        expect(other_role_input).to_be_visible()
        other_role_input.fill(new_role)

    page.get_by_role("button", name="Update Profile").click()
    expect(page).to_have_url(f"{FRONTEND_BASE_URL}/profile/")

    expect(page.get_by_text(new_name)).to_be_visible()
    expect(page.get_by_text(new_organisation)).to_be_visible()
    expect(page.get_by_text(new_role)).to_be_visible()


@pytest.mark.parametrize(
    "new_name, new_organisation, new_role",
    [
        ("Updated Profile User", "Uppsala University", "Postdoctoral Researcher"),
        ("Another User", "Other Organisation", "Other Role"),
        ("Testing User", "Lund University", "PhD Student"),
        ("Onenameuser", "Lund University", "Non-dropdown Role"),
        ("Test Test User", "Other Organisation", "PhD Student"),
    ],
)
def test_edit_profile(
    logged_in_edit_profile_user_page: Page,
    new_name: str,
    new_organisation: str,
    new_role: str,
):
    """
    Test editing the profile of a user, including handling of "Other" option for organisation and role.
    """
    _edit_profile(
        page=logged_in_edit_profile_user_page,
        new_name=new_name,
        new_organisation=new_organisation,
        new_role=new_role,
    )


def test_edit_profile_cancel_does_not_save_changes(logged_in_edit_profile_user_page: Page):
    """
    Test canceling the profile edit with changes made on the form does not save the changes.
    """
    page = logged_in_edit_profile_user_page

    navigate_to(page=page, path="/profile/edit")

    page.get_by_role("textbox", name="Full Name").fill("Cancelled name")
    page.get_by_label("Organisation", exact=True).select_option("Other")
    page.get_by_role("textbox", name="Your organisation").fill("Cancelled Organisation")
    page.get_by_label("Role", exact=True).select_option("Other")
    page.get_by_role("textbox", name="Your Role").fill("Cancelled Role")

    page.get_by_role("link", name="Cancel").click()
    expect(page).to_have_url(f"{FRONTEND_BASE_URL}/profile/")

    for cancelled_text in ["Cancelled Update", "Cancelled Organisation", "Cancelled Role"]:
        expect(page.get_by_text(cancelled_text)).not_to_be_visible()
