"""
unit tests for displaying announcements
"""

from unittest.mock import MagicMock, call, patch

from divbase_cli.services.announcements import get_and_display_announcements


@patch("divbase_cli.services.announcements.make_unauthenticated_request")
@patch("divbase_cli.services.announcements.print")
def test_get_and_display_announcements_no_announcements(mock_print: MagicMock, mock_request: MagicMock) -> None:
    """Test that nothing is printed if the API returns an empty list."""
    mock_response = MagicMock()
    mock_response.json.return_value = []
    mock_request.return_value = mock_response

    get_and_display_announcements("https://divbase.com")

    mock_request.assert_called_once_with(
        method="GET",
        divbase_base_url="https://divbase.com",
        api_route="/v1/core/announcements",
    )
    mock_print.assert_not_called()


@patch("divbase_cli.services.announcements.make_unauthenticated_request")
@patch("divbase_cli.services.announcements.print")
def test_get_and_display_announcements_single_announcement(mock_print: MagicMock, mock_request: MagicMock) -> None:
    """Test that a single announcement is printed correctly."""
    mock_response = MagicMock()
    mock_response.json.return_value = [
        {"heading": "Test Heading", "message": "Test Message", "level": "info"},
    ]
    mock_request.return_value = mock_response

    get_and_display_announcements("https://divbase.com")

    expected_calls = [
        call("[bold blue]Test Heading[/bold blue]"),
        call("Test Message\n"),
    ]
    mock_print.assert_has_calls(expected_calls)


@patch("divbase_cli.services.announcements.make_unauthenticated_request")
@patch("divbase_cli.services.announcements.print")
def test_get_and_display_announcements_multiple_announcements(mock_print: MagicMock, mock_request: MagicMock) -> None:
    """Test that multiple announcements are printed correctly as a numbered list."""
    mock_response = MagicMock()
    mock_response.json.return_value = [
        {"heading": "Heading 1", "message": "Message 1", "level": "success"},
        {"heading": "Heading 2", "message": "Message 2", "level": "danger"},
    ]
    mock_request.return_value = mock_response

    get_and_display_announcements("https://divbase.com")

    expected_calls = [
        call("[bold green]1) Heading 1[/bold green]"),
        call("Message 1\n"),
        call("[bold red]2) Heading 2[/bold red]"),
        call("Message 2\n"),
    ]
    mock_print.assert_has_calls(expected_calls)


@patch("divbase_cli.services.announcements.make_unauthenticated_request")
@patch("divbase_cli.services.announcements.print")
def test_get_and_display_announcements_unexpected_level(mock_print: MagicMock, mock_request: MagicMock) -> None:
    """Test that an announcement with an unexpected level defaults to no color."""
    mock_response = MagicMock()
    mock_response.json.return_value = [
        {"heading": "Weird Level", "message": "No color message", "level": "critical"},
    ]
    mock_request.return_value = mock_response

    get_and_display_announcements("https://divbase.com")

    # The color mapping will return None, resulting in rich using default styling.
    expected_calls = [
        call("[bold None]Weird Level[/bold None]"),
        call("No color message\n"),
    ]
    mock_print.assert_has_calls(expected_calls)
