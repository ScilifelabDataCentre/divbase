"""
Test register page functionality using Playwright.
"""

import re

from playwright.sync_api import Page, expect

from .conftest import FRONTEND_BASE_URL, is_logged_in_as, login_via_login_form, navigate_to

GENERIC_REGISTRATION_ERROR_MESSAGE = "Registration failed, please try again."
EMAIL_NOT_VERIFIED_MESSAGE = "Email address not verified"


def test_register_new_user_complete_flow(page: Page, mailpit_page: Page):
    """Test complete user registration flow including email verification step."""
    test_name = "New Test User"
    test_email = "newuser@gmail.com"
    test_password = "badpassword"

    navigate_to(page, "/register")
    expect(page).to_have_url(f"{FRONTEND_BASE_URL}/register")

    page.get_by_role("textbox", name="Full Name").fill(test_name)
    page.get_by_role("textbox", name="Email Address").fill(test_email)
    page.get_by_role("textbox", name="Password", exact=True).fill(test_password)
    page.get_by_role("textbox", name="Confirm Password").fill(test_password)
    page.get_by_role("button", name=" Create Account").click()

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
    navigate_to(page, "/register")

    existing_email = EXISTING_ACCOUNTS["ADMIN_USER"]["email"]

    page.get_by_role("textbox", name="Full Name").fill("Test User")
    page.get_by_role("textbox", name="Email Address").fill(existing_email)
    page.get_by_role("textbox", name="Password", exact=True).fill("newpassword")
    page.get_by_role("textbox", name="Confirm Password").fill("newpassword")
    page.get_by_role("button", name=" Create Account").click()

    expect(page.get_by_text(GENERIC_REGISTRATION_ERROR_MESSAGE)).to_be_visible()


def test_register_navigation_links(page: Page):
    """Test navigation links on register page."""
    navigate_to(page, "/register")

    # Target main to avoid navbar links for these pages.
    login_button = page.get_by_role("main").get_by_role("link", name=re.compile("Login"))
    expect(login_button).to_be_visible()

    forgot_password_button = page.get_by_role("main").get_by_role("link", name=re.compile("Reset Password"))
    expect(forgot_password_button).to_be_visible()
