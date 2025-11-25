""" """

from playwright.sync_api import Page, expect

from .conftest import FRONTEND_BASE_URL, is_logged_in_as, login_via_login_form, navigate_to


def test_admin_login_success(CREDENTIALS, page: Page):
    """Test successful login with admin credentials."""
    navigate_to(page, "/login")
    expect(page).to_have_url(f"{FRONTEND_BASE_URL}/login")

    expect(page.locator('input[type="email"], input[name="email"]')).to_be_visible()
    expect(page.locator('input[type="password"], input[name="password"]')).to_be_visible()

    admin_email = CREDENTIALS["ADMIN_USER"]["email"]
    admin_password = CREDENTIALS["ADMIN_USER"]["password"]
    login_via_login_form(page=page, email=admin_email, password=admin_password)

    assert is_logged_in_as(page=page, email=admin_email)


def test_edit_user_login_success(CREDENTIALS, page: Page):
    """Test successful login with edit user credentials."""
    navigate_to(page, "/login")
    expect(page).to_have_url(f"{FRONTEND_BASE_URL}/login")

    expect(page.locator('input[type="email"], input[name="email"]')).to_be_visible()
    expect(page.locator('input[type="password"], input[name="password"]')).to_be_visible()

    edit_user_email = CREDENTIALS["EDIT_USER"]["email"]
    edit_user_password = CREDENTIALS["EDIT_USER"]["password"]
    login_via_login_form(page=page, email=edit_user_email, password=edit_user_password)

    assert is_logged_in_as(page=page, email=edit_user_email)
