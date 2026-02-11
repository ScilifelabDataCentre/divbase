"""
Test register page functionality using Playwright.
"""

import re

import pytest
from playwright.sync_api import Page, expect

from .conftest import is_logged_in_as, login_via_login_form, navigate_to, register_new_user

GENERIC_REGISTRATION_ERROR_MESSAGE = "Registration failed, please try again."
EMAIL_NOT_VERIFIED_MESSAGE = "Email address not verified"


@pytest.mark.parametrize(
    "test_name,test_email,test_organisation,test_role",
    [
        (
            "New Test User Uppsala",
            "newuser@uu.se",
            "Uppsala University",
            "Postdoctoral Researcher",
        ),
        (
            "New Test User EMBL",
            "newuser@embl.de",
            "European Bioinformatics Institute (EMBL-EBI)",
            "Developer",
        ),
    ],
)
def test_register_new_user_complete_flow(
    page: Page,
    mailpit_page: Page,
    test_name: str,
    test_email: str,
    test_organisation: str,
    test_role: str,
):
    """
    Test complete user registration flow including email verification step.

    We test creating a user from a Swedish university (in the organisation dropdown)
    and a user from an organisation not in the dropdown to ensure both flows work correctly.
    """
    test_password = "badpassword"
    register_new_user(
        page=page,
        name=test_name,
        email=test_email,
        organisation=test_organisation,
        role=test_role,
        password=test_password,
    )

    # validate logging in now does not work
    navigate_to(page, "/login")
    login_via_login_form(page, test_email, test_password)
    expect(page.get_by_text(EMAIL_NOT_VERIFIED_MESSAGE)).to_be_visible()

    # Check that verification email is received, open it and click on verification link
    email_link = mailpit_page.get_by_role("link", name=f"DivBase To: {test_email}")
    expect(email_link).to_be_visible()
    email_link.click()

    verification_link = mailpit_page.locator("#preview-html").content_frame.get_by_role(
        "link", name="Verify Email Address"
    )

    with mailpit_page.expect_popup() as new_tab_info:
        verification_link.click()

    # Back on divbase site, click verify email and test can login.
    new_tab = new_tab_info.value
    try:
        new_tab.get_by_role("button", name="Confirm My Email Address").click()
        login_via_login_form(new_tab, test_email, test_password)
        assert is_logged_in_as(new_tab, test_email)
    finally:
        new_tab.close()


def test_register_with_existing_email(page: Page, EXISTING_ACCOUNTS):
    """Test registration with email that already exists, does not tell user that the email already exists."""
    existing_email = EXISTING_ACCOUNTS["ADMIN_USER"]["email"]
    register_new_user(
        page=page,
        name="Test User",
        email=existing_email,
        organisation="Uppsala University",
        role="Lecturer",
        password="newpassword",
        expect_success=False,
    )
    expect(page.get_by_text(GENERIC_REGISTRATION_ERROR_MESSAGE)).to_be_visible()


def test_register_navigation_links(page: Page):
    """Test navigation links on register page."""
    navigate_to(page, "/register")

    # Target main to avoid navbar links for these pages.
    login_button = page.get_by_role("main").get_by_role("link", name=re.compile("Login"))
    expect(login_button).to_be_visible()

    forgot_password_button = page.get_by_role("main").get_by_role("link", name=re.compile("Reset Password"))
    expect(forgot_password_button).to_be_visible()
