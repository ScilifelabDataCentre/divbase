"""
Send automatic emails for e.g. email verification, password reset, notifications, etc.

In production/when deployed we rely on using smtp-relay.gmail.com

For local development/testing we have a Mailpit service that is ran as part of the docker compose stack.
All email sent by the API will be caught by Mailpit and can be viewed in the Mailpit web UI at http://localhost:8025

### Templates:
Email templates are written using MJML template system and converted to HTML.
The compiled HTML templates are used by Jinja2 to create the emails actual content.
To compile the templates I used the vscode "MJML Official" extension.
"""

import logging
from pathlib import Path
from time import sleep
from typing import Any

import emails
from jinja2 import Template

from divbase_api.api_config import api_settings
from divbase_api.security import TokenType, create_token

logger = logging.getLogger(__name__)


def render_email_template(template_name: str, context: dict[str, Any]) -> str:
    """
    Render an email template with Jinja2.
    (Jinja2 relies on the context dict to fill in the variables in the template.)
    """
    context["support_email"] = api_settings.general.user_support_email
    email_templates = Path(__file__).parent / "email_templates" / "build"
    template_str = (email_templates / template_name).read_text()
    return Template(template_str).render(context)


def _send_email(email_to: str, subject: str, html_content: str, retries: int = 3, retry_delay: int = 10) -> None:
    """
    Helper function to send any type of email,
    use one of the specific email functions below instead.
    """
    message = emails.Message(
        subject=subject,
        html=html_content,
        mail_from=("DivBase", api_settings.email.from_email),
    )
    smtp_options = {"host": api_settings.email.smtp_server, "port": api_settings.email.smtp_port}

    if api_settings.email.smtp_tls:
        smtp_options["tls"] = True
    elif api_settings.email.smtp_ssl:
        smtp_options["ssl"] = True

    if api_settings.email.smtp_user:
        smtp_options["user"] = api_settings.email.smtp_user
    if api_settings.email.smtp_password:
        smtp_options["password"] = api_settings.email.smtp_password.get_secret_value()

    for attempt in range(1, retries + 1):
        response = message.send(to=email_to, smtp=smtp_options)

        if response.status_code == 250:
            logger.info(f"Email to {email_to} with subject '{subject}' sent successfully.")
            return
        else:
            logger.warning(f"Failed to send email to {email_to} with subject '{subject}'. Response: {response}")

        if attempt < retries:
            sleep(retry_delay)

    logger.error(f"All {retries} attempts to send email to {email_to} with subject '{subject}', failed.")


def send_test_email(email_to: str) -> None:
    """
    Send a test email to the specified email address.
    """
    subject = "DivBase - test email"
    html_content = render_email_template(
        template_name="test_email.html",
        context={"email": email_to},
    )
    _send_email(email_to=email_to, subject=subject, html_content=html_content)


def send_verification_email(email_to: str, user_id: int) -> None:
    """
    Send a verification email to the specified email address.
    """
    subject = "DivBase - verify your email address"

    token_data = create_token(subject=user_id, token_type=TokenType.EMAIL_VERIFICATION)
    verification_url = f"{api_settings.general.frontend_base_url}/verify-email?token={token_data.token}"

    link_expire_hours = api_settings.email.email_verify_expires_seconds // 3600

    html_content = render_email_template(
        template_name="email_verification.html",
        context={"email": email_to, "verification_url": verification_url, "link_expire_hours": link_expire_hours},
    )
    _send_email(email_to=email_to, subject=subject, html_content=html_content)


def send_email_already_verified_email(email_to: str) -> None:
    """
    Send an email informing the user that their email is already verified.

    Hit if user requests a new verification email but their email is already verified.
    """
    subject = "DivBase - Your email is already verified"
    html_content = render_email_template(
        template_name="email_already_verified.html",
        context={"email": email_to},
    )
    _send_email(email_to=email_to, subject=subject, html_content=html_content)


def send_password_reset_email(email_to: str, user_id: int) -> None:
    """
    Send a password reset email to the specified email address.
    """
    subject = "DivBase - reset your password"

    token_data = create_token(subject=user_id, token_type=TokenType.PASSWORD_RESET)
    reset_password_url = f"{api_settings.general.frontend_base_url}/reset-password?token={token_data.token}"

    link_expire_hours = api_settings.email.password_reset_expires_seconds // 3600

    html_content = render_email_template(
        template_name="reset_password.html",
        context={"email": email_to, "password_reset_url": reset_password_url, "link_expire_hours": link_expire_hours},
    )
    _send_email(email_to=email_to, subject=subject, html_content=html_content)


def send_password_has_been_reset_email(email_to: str) -> None:
    """
    Send a email to tell the user their password has now been reset.
    """
    subject = "DivBase - Your password has been reset"

    html_content = render_email_template(
        template_name="password_was_reset.html",
        context={"email": email_to},
    )
    _send_email(email_to=email_to, subject=subject, html_content=html_content)


def send_pat_created_email(email_to: str, pat_name: str, expires_at_cet: str | None) -> None:
    """
    Send a security notification email when a new personal access token is created.
    """
    subject = "DivBase - New personal access token created"
    pats_url = f"{api_settings.general.frontend_base_url}/pats"
    html_content = render_email_template(
        template_name="pat_created.html",
        context={"pat_name": pat_name, "pats_url": pats_url, "expires_at_cet": expires_at_cet},
    )
    _send_email(email_to=email_to, subject=subject, html_content=html_content)


def send_pat_revoked_email(email_to: str, pat_name: str) -> None:
    """
    Send a security notification email when a personal access token is revoked.
    """
    subject = "DivBase - Personal access token revoked"
    pats_url = f"{api_settings.general.frontend_base_url}/pats"
    html_content = render_email_template(
        template_name="pat_revoked.html",
        context={"pat_name": pat_name, "pats_url": pats_url},
    )
    _send_email(email_to=email_to, subject=subject, html_content=html_content)
