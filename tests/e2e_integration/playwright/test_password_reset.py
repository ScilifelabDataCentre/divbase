"""
Tests for the password reset functionality.

Users can reset password via email link.
They can get a reset link when logged in by going to their profile page, OR
when not logged in via the login/register pages (buttons for forgot password).
"""

import re

from playwright.sync_api import Page, expect

from divbase_api.frontend_routes.auth import INVALID_EXPIRED_PASSWORD_TOKEN_MSG

from .conftest import (
    FRONTEND_BASE_URL,
    MAILPIT_BASE_URL,
    is_logged_in_as,
    login_via_login_form,
    navigate_to,
)
from .test_login_page import LOGIN_FAILED_MESSAGE

PASSWORD_RESET_LINK_SENT_ALERT_MESSAGE = re.compile(
    "If your account exists, and your email is verified, a password reset email has been sent to *"
)
PASSWORD_RESET_ALERT_MESSAGE = "Your password has been reset successfully, you can now log in."


def reset_password_flow(page: Page, mailpit_page: Page, user_email: str, new_password: str):
    """Complete password reset flow steps starting from the reset password page"""
    page.get_by_placeholder("Enter your email address").fill(user_email)
    page.get_by_role("button", name=re.compile("Send Reset Link")).click()
    reset_alert = page.get_by_role("alert")
    expect(reset_alert).to_have_text(PASSWORD_RESET_LINK_SENT_ALERT_MESSAGE)

    # Check that reset email is received, open it and click on reset link
    email_link = mailpit_page.get_by_role("link", name=f"To: {user_email} DivBase - reset your password")
    expect(email_link).to_be_visible()
    email_link.click()

    reset_link = mailpit_page.locator("#preview-html").content_frame.get_by_role("link", name="Reset Password")
    with mailpit_page.expect_popup() as new_tab_info:
        reset_link.click()
    new_tab = new_tab_info.value

    # reset password in the new tab
    try:
        new_tab.get_by_placeholder("Enter your new password").fill(new_password)
        new_tab.get_by_placeholder("Confirm your new password").fill(new_password)
        new_tab.get_by_role("button", name="Reset Password").click()
        reset_alert = new_tab.get_by_role("alert")
        expect(reset_alert).to_have_text(PASSWORD_RESET_ALERT_MESSAGE)
    finally:
        new_tab.close()

    # validate we recieved a "your password has been reset" email
    password_has_been_reset_link = mailpit_page.get_by_role(
        "link", name=f"To: {user_email} DivBase - Your password has been reset"
    )
    expect(password_has_been_reset_link).to_be_visible()


def test_password_reset_flow(page: Page, mailpit_page: Page, EXISTING_ACCOUNTS):
    """Test complete password reset flow for a user."""
    user_email = EXISTING_ACCOUNTS["FORGOT_PASSWORD_USER"]["email"]
    current_password = EXISTING_ACCOUNTS["FORGOT_PASSWORD_USER"]["password"]
    new_password = "badpassword123"

    navigate_to(page, "/login")
    page.get_by_role("link", name=re.compile("Reset Password")).click()

    reset_password_flow(
        page=page,
        mailpit_page=mailpit_page,
        user_email=user_email,
        new_password=new_password,
    )

    # validate can't login with old password but can with new password
    navigate_to(page, "/login")
    login_via_login_form(page, user_email, current_password)
    expect(page.get_by_text(LOGIN_FAILED_MESSAGE)).to_be_visible()
    expect(page).to_have_url(f"{FRONTEND_BASE_URL}/login")

    login_via_login_form(page, user_email, new_password)
    assert is_logged_in_as(page=page, email=user_email)


def test_password_reset_for_non_existent_user(page: Page, mailpit_page: Page):
    """
    The site should not reveal whether an email exists in the system or not.
    but no email should be sent and should not be clear to user that the email does not exist in the system.
    """
    navigate_to(page, "/login")
    false_user_email = "non-existent-user@example.com"

    page.get_by_role("link", name=re.compile("Reset Password")).click()
    page.get_by_placeholder("Enter your email address").fill(false_user_email)
    page.get_by_role("button", name=re.compile("Send Reset Link")).click()
    reset_alert = page.get_by_role("alert")
    expect(reset_alert).to_have_text(PASSWORD_RESET_LINK_SENT_ALERT_MESSAGE)

    email_link = mailpit_page.get_by_role("link", name=re.compile("DivBase - reset your password"))
    expect(email_link).to_have_count(0)


def test_password_reset_as_logged_in_user(page: Page, mailpit_page: Page, EXISTING_ACCOUNTS):
    """
    After password reset, user should be logged out and needs to login with new password.
    User accesses the password reset flow while logged in via their profile page.

    Need to validate that after password reset, user is logged out and needs to login with new password.
    """
    user_email = EXISTING_ACCOUNTS["PASSWORD_RESET_USER"]["email"]
    current_password = EXISTING_ACCOUNTS["PASSWORD_RESET_USER"]["password"]
    new_password = "badpassword123"

    login_via_login_form(page=page, email=user_email, password=current_password)

    navigate_to(page, "/profile")
    page.get_by_role("link", name="Request Password Reset Email").click()

    reset_password_flow(
        page=page,
        mailpit_page=mailpit_page,
        user_email=user_email,
        new_password=new_password,
    )
    navigate_to(page, "/")
    page.reload()
    expect(page.get_by_role("link", name=re.compile("Log in"))).to_be_visible()


def test_password_reset_token_cannot_be_reused(page: Page, mailpit_page: Page, EXISTING_ACCOUNTS):
    """
    After password reset, should not be able to reuse the same token again.
    Due to token being revoked after use.
    """
    user_email = EXISTING_ACCOUNTS["FORGOT_PASSWORD_USER"]["email"]
    new_password = "badpassword123"

    navigate_to(page, "/login")
    page.get_by_role("link", name=re.compile("Reset Password")).click()

    reset_password_flow(
        page=page,
        mailpit_page=mailpit_page,
        user_email=user_email,
        new_password=new_password,
    )

    # Try to reuse the same reset link again
    mailpit_page.goto(MAILPIT_BASE_URL)
    email_link = mailpit_page.get_by_role("link", name=f"To: {user_email} DivBase - reset your password")
    expect(email_link).to_be_visible()
    email_link.click()

    reset_link = mailpit_page.locator("#preview-html").content_frame.get_by_role("link", name="Reset Password")
    with mailpit_page.expect_popup() as new_tab_info:
        reset_link.click()
    new_tab = new_tab_info.value

    # validate we see error about invalid/expired token and have option to request a new link.
    expect(new_tab.get_by_text(INVALID_EXPIRED_PASSWORD_TOKEN_MSG)).to_be_visible()
    expect(new_tab.get_by_placeholder("Enter your email address")).to_be_visible()
