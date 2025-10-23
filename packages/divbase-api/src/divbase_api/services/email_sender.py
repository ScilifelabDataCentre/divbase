"""
Send email via an SMTP server.

In production/deployed we rely on TODO - explain when figured out...

For local development/testing we have an optional Mailpit service that can be run as part of the docker compose stack.
If running, all emails will be caught by Mailpit and can be viewed in the Mailpit web UI at http://localhost:8025

### Templates:
Email templates are written using MJML template system and converted to HTML using the MJML CLI tool.
The compiled HTML templates are used by Jinja2 to create the emails actual content.
To compile the templates I used the vscode "MJML Official" extension.
"""

import logging
from pathlib import Path
from typing import Any

import emails
from jinja2 import Template

from divbase_api.config import settings

logger = logging.getLogger(__name__)


EMAIL_TEMPLATES = Path(__file__).parent / "email_templates" / "build"


def render_email_template(template_name: str, context: dict[str, Any]) -> str:
    """
    Render an email template with Jinja2.
    (Jinja2 relies on the context dict to fill in the variables in the template.)
    """
    template_str = (EMAIL_TEMPLATES / template_name).read_text()
    return Template(template_str).render(context)


def _send_email(email_to: str, subject: str, html_content: str) -> None:
    """
    Helper function to send any type of email,
    use one of the specific email functions below instead.
    TODO: Some kind of guard clause if emails are disabled in settings?
    """
    message = emails.Message(
        subject=subject,
        html=html_content,
        mail_from=("DivBase", settings.email.from_email),
    )
    smtp_options = {"host": settings.email.smtp_server, "port": settings.email.smtp_port}

    if settings.email.smtp_tls:
        smtp_options["tls"] = True
    elif settings.email.smtp_ssl:
        smtp_options["ssl"] = True

    if settings.email.smtp_user:
        smtp_options["user"] = settings.email.smtp_user
    if settings.email.smtp_password:
        smtp_options["password"] = settings.email.smtp_password.get_secret_value()

    response = message.send(to=email_to, smtp=smtp_options)
    logger.info(f"Email sent with response: {response}")


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


# TODO
def send_verification_email():
    pass


def send_password_reset_email():
    pass
