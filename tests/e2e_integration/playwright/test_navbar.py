"""
Tests for the navbar component using Playwright

Tests run on home page as navbar is same on all pages.
"""

from playwright.sync_api import Page, expect

from .conftest import FRONTEND_BASE_URL, navigate_to


def test_admin_can_access_admin_panel(logged_in_admin_page: Page):
    """Test admin user can access admin panel."""
    navigate_to(logged_in_admin_page, "/admin")
    expect(logged_in_admin_page).to_have_url(f"{FRONTEND_BASE_URL}/admin/")


def test_non_admin_cannot_access_admin_panel(
    logged_in_manager_page: Page, logged_in_edit_page: Page, logged_in_read_page: Page
):
    """Test non-admin users types cannot access admin panel."""
    navigate_to(logged_in_manager_page, "/admin")
    expect(logged_in_manager_page).to_have_url(f"{FRONTEND_BASE_URL}/")

    navigate_to(logged_in_edit_page, "/admin")
    expect(logged_in_edit_page).to_have_url(f"{FRONTEND_BASE_URL}/")

    navigate_to(logged_in_read_page, "/admin")
    expect(logged_in_read_page).to_have_url(f"{FRONTEND_BASE_URL}/")


def test_non_logged_in_user_cannot_access_admin_panel(page: Page):
    """Test navbar for non-logged in user."""
    navigate_to(page, "/admin")
    expect(page).to_have_url(f"{FRONTEND_BASE_URL}/")
