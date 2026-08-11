"""
Functions that resolve for the CLI commands things like:
    - which project to use
    - which download directory to use
    - which DivBase API URL to use
Based on provided user input and their config file.
"""

from pathlib import Path

from divbase_cli.cli_config import cli_settings
from divbase_cli.cli_exceptions import AuthenticationError, ProjectNameNotSpecifiedError
from divbase_cli.user_auth import get_pat_for_authentication
from divbase_cli.user_config import ProjectConfig, load_user_config


def resolve_and_authenticate_project(project_name: str | None) -> ProjectConfig:
    """
    Helper function that resolves the project to use for a CLI command and ensures
    the user is logged in to the correct DivBase API URL for that project.

    Returns the ProjectConfig object for the resolved project, will raise errors otherwise.
    """
    user_config = load_user_config()
    if not project_name:
        project_name = user_config.default_project
    if not project_name:
        raise ProjectNameNotSpecifiedError()

    project_config = user_config.project_info(project_name)

    if user_config.logged_in_url:
        if project_config.divbase_url != user_config.logged_in_url:
            raise AuthenticationError(
                f"You are trying to run a command for the project '{project_name}' which has the DivBase URL: {project_config.divbase_url} \n"
                f"You are however logged in to this URL: {user_config.logged_in_url} \n"
                "Please log in again."
            )
        return project_config

    if get_pat_for_authentication():
        # if a user has a PAT, we will assume the PAT is for the correct URL
        return project_config

    raise AuthenticationError("You are not logged in. Please log in with 'divbase-cli auth login [EMAIL]'.")


def resolve_url_for_non_project_specific_commands() -> str:
    """
    Resolve the DivBase API URL to use for CLI commands that are not project-specific.

    Returns the url the user is either logged into or in the case of a user using a personal access token (PAT),
    the default API URL, since PATs don't require login to be used.
    Current examples: auth whoami and some task-history commands that are not project-specific.

    Priority mirrors make_authenticated_request: active session first, PAT as fallback.
    """
    config = load_user_config()
    if config.logged_in_url:
        return config.logged_in_url

    if get_pat_for_authentication():
        return cli_settings.DIVBASE_API_URL

    raise AuthenticationError("You are not logged in. Please log in with 'divbase-cli auth login [EMAIL]'.")


def resolve_download_dir(download_dir: str | None) -> Path:
    """
    Helper function to resolve the download directory to use for a CLI command involving downloading files.

    Priority given to `download_dir` argument, then if a default is set in the user config.
    Note: "." or None should default to the current working directory.
    """
    if not download_dir:
        config = load_user_config()
        download_dir = config.default_download_dir

    if download_dir and download_dir != ".":
        return Path(download_dir).expanduser()
    return Path.cwd()
