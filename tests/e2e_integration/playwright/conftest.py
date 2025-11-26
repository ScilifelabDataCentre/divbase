"""
Fixtures and helper functions for Playwright tests, for testing the DivBase frontend
"""

import re

import pytest
from playwright.sync_api import BrowserContext, Page

FRONTEND_BASE_URL = "http://localhost:8001"
MAILPIT_BASE_URL = "http://localhost:8026"


def navigate_to(page: Page, path: str):
    """Navigate to a specific page"""
    page.goto(f"{FRONTEND_BASE_URL}{path}")


def login_via_login_form(page: Page, email: str, password: str):
    page.get_by_role("textbox", name="Email Address").fill(email)
    page.get_by_role("textbox", name="Password").fill(password)
    page.get_by_role("button", name=re.compile(r"Sign In", re.IGNORECASE)).click()


def is_logged_in_as(page: Page, email: str) -> bool:
    """Check if the user is logged in as the specified email."""
    return page.get_by_role("button", name=f"User menu for {email}").count() > 0


@pytest.fixture(scope="session")
def EXISTING_ACCOUNTS(CONSTANTS):
    """
    Credentials for tests.
    Just reformats those in CONSTANTS for easier access here as will be used so much here.
    """
    return {
        "ADMIN_USER": CONSTANTS["ADMIN_CREDENTIALS"],
        "READ_USER": CONSTANTS["TEST_USERS"]["read user"],
        "EDIT_USER": CONSTANTS["TEST_USERS"]["edit user"],
        "MANAGE_USER": CONSTANTS["TEST_USERS"]["manage user"],
    }


@pytest.fixture
def mailpit_page(context: BrowserContext):
    """Fixture to provide a Mailpit page that auto-closes after each test."""
    page = context.new_page()
    page.goto(MAILPIT_BASE_URL)
    yield page

    # delete all emails after test to keep Mailpit clean
    page.get_by_role("button", name="Delete all").click()
    page.close()


def _create_logged_in_user_page(context: BrowserContext, user_credentials: dict) -> Page:
    """Helper to create a logged-in user page."""
    page = context.new_page()
    navigate_to(page, "/login")
    login_via_login_form(page, user_credentials["email"], user_credentials["password"])

    assert is_logged_in_as(page, user_credentials["email"]), f"Failed to log in as {user_credentials['email']}"

    page.goto(FRONTEND_BASE_URL)
    return page


@pytest.fixture
def logged_in_admin_page(context: BrowserContext, EXISTING_ACCOUNTS):
    """Logged in admin user at home page."""
    page = _create_logged_in_user_page(context, EXISTING_ACCOUNTS["ADMIN_USER"])
    yield page
    page.close()


@pytest.fixture
def logged_in_manager_page(context: BrowserContext, EXISTING_ACCOUNTS):
    """Logged in manager user at home page."""
    page = _create_logged_in_user_page(context, EXISTING_ACCOUNTS["MANAGE_USER"])
    yield page
    page.close()


@pytest.fixture
def logged_in_edit_page(context: BrowserContext, EXISTING_ACCOUNTS):
    """Logged in edit user at home page."""
    page = _create_logged_in_user_page(context, EXISTING_ACCOUNTS["EDIT_USER"])
    yield page
    page.close()


@pytest.fixture
def logged_in_read_page(context: BrowserContext, EXISTING_ACCOUNTS):
    """Logged in read user at home page."""
    page = _create_logged_in_user_page(context, EXISTING_ACCOUNTS["READ_USER"])
    yield page
    page.close()
