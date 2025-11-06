"""
Unit tests for the email sender service.
"""

from unittest.mock import ANY, patch

from divbase_api.services.email_sender import (
    _send_email,
    render_email_template,
    send_test_email,
    send_verification_email,
)


@patch("emails.Message.send")
def test_send_email_success(mock_send):
    mock_send.return_value.status_code = 250

    _send_email(
        email_to="test@example.com",
        subject="Test Subject",
        html_content="<p>Test Content</p>",
    )

    mock_send.assert_called_once()


@patch("divbase_api.services.email_sender._send_email")
def test_send_test_email(mock_send_email):
    send_test_email(email_to="test@example.com")

    expected_html_content = render_email_template(
        template_name="test_email.html",
        context={"email": "test@example.com"},
    )

    mock_send_email.assert_called_once_with(
        email_to="test@example.com",
        subject="DivBase - test email",
        html_content=expected_html_content,
    )


@patch("divbase_api.services.email_sender._send_email")
def test_send_verification_email(mock_send_email):
    send_verification_email(email_to="test@example.com", user_id=123)

    mock_send_email.assert_called_once_with(
        email_to="test@example.com",
        subject="DivBase - verify your email address",
        html_content=ANY,  # too large to compare directly
    )
