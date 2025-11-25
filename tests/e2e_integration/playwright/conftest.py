"""
Fixtures and helper functions for Playwright tests, for testing the DivBase frontend
"""

import re

import pytest
from playwright.sync_api import Page

FRONTEND_BASE_URL = "http://localhost:8001"


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
def CREDENTIALS(CONSTANTS):
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
