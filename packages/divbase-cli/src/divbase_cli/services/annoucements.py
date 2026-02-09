"""
Get announcements from divbase API and display them to the user.
"""

from rich import print

from divbase_cli.user_auth import make_authenticated_request
from divbase_lib.api_schemas.annoucements import AnnouncementResponse

COLOR_MAPPING = {
    "info": "blue",
    "success": "green",
    "warning": "yellow",
    "danger": "red",
}


def get_and_display_announcements(divbase_base_url: str) -> None:
    """
    Get announcements from the API and display them to the user.
    Used when logging in to display any important messages to the user.
    """
    response = make_authenticated_request(
        method="GET",
        divbase_base_url=divbase_base_url,
        api_route="/v1/core/announcements",
    )
    announcements_data = response.json()
    if not announcements_data:
        return
    announcements = [AnnouncementResponse.model_validate(ann) for ann in announcements_data]

    print("[bold]Announcements from DivBase Server:[/bold]")

    if len(announcements) == 1:
        color = COLOR_MAPPING[announcements[0].level]
        print(f"[bold {color}]{announcements[0].heading}[/bold {color}]")
        print(f"{announcements[0].message}\n")
        return

    for idx, ann in enumerate(announcements, start=1):
        color = COLOR_MAPPING[ann.level]
        print(f"[bold {color}]{idx}) {ann.heading}[/bold {color}]")
        print(f"{ann.message}\n")
