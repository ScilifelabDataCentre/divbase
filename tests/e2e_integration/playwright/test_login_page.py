""" """

import re

import pytest
from playwright.sync_api import Page, expect

from .conftest import FRONTEND_BASE_URL, is_logged_in_as, login_via_login_form, navigate_to

LOGIN_FAILED_MESSAGE = "Invalid email or password or user account does not exist"


@pytest.mark.parametrize("user_type", ["ADMIN_USER", "READ_USER", "EDIT_USER", "MANAGE_USER"])
def test_user_login_success(CREDENTIALS, page: Page, user_type: str):
    """Test successful login for all user types."""
    navigate_to(page, "/login")
    expect(page).to_have_url(f"{FRONTEND_BASE_URL}/login")

    expect(page.locator('input[type="email"], input[name="email"]')).to_be_visible()
    expect(page.locator('input[type="password"], input[name="password"]')).to_be_visible()

    user_email = CREDENTIALS[user_type]["email"]
    user_password = CREDENTIALS[user_type]["password"]

    login_via_login_form(page=page, email=user_email, password=user_password)
    assert is_logged_in_as(page=page, email=user_email)


def test_invalid_login_credentials(page: Page):
    """Test login with invalid credentials shows error."""
    navigate_to(page, "/login")
    expect(page.get_by_text(LOGIN_FAILED_MESSAGE)).not_to_be_visible()
    login_via_login_form(page=page, email="invalid@email.com", password="wrongpassword")
    expect(page.get_by_text(LOGIN_FAILED_MESSAGE)).to_be_visible()
    expect(page).to_have_url(f"{FRONTEND_BASE_URL}/login")


def test_navigation_links_present(page: Page):
    """Test that navigation links to other auth pages are present."""
    navigate_to(page, "/login")

    # Target main to avoid navbar links for these pages.
    register_button = page.get_by_role("main").get_by_role("link", name=re.compile("Create Account"))
    expect(register_button).to_be_visible()

    forgot_password_button = page.get_by_role("main").get_by_role("link", name=re.compile("Reset Password"))
    expect(forgot_password_button).to_be_visible()


# TODO - important tests to do, but require thought about setup.
# test_login with user with disabled account
# test_login with unverified email - but wrong password
# test_login with unverified email - correct password - gets better message
