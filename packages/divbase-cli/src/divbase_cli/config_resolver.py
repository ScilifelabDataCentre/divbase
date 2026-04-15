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
from divbase_cli.user_config import ProjectConfig, load_user_config


def ensure_logged_in(desired_url: str | None = None) -> str:
    """
    Ensure the user is logged in by checking the logged_in_url value in the user config.

    Optionally checks the logged_in_url matches the desired_url (e.g. for project-specific commands) if provided.
    """
    # when using a personal access token (PAT), we can skip the login check
    # and use either project's URL.
    # For commands that don't use a or the default API URL
    if cli_settings.DIVBASE_API_PAT:
        return desired_url or cli_settings.DIVBASE_API_URL

    config = load_user_config()
    if not config.logged_in_url:
        raise AuthenticationError("You are not logged in. Please log in with 'divbase-cli auth login [EMAIL]'.")
    if desired_url and config.logged_in_url != desired_url:
        raise AuthenticationError(
            f"You are not logged in to the correct DivBase URL: {desired_url}. Please log in again."
        )
    return config.logged_in_url


def resolve_url_for_non_project_specific_commands() -> str:
    """
    Resolve the DivBase API URL to use for CLI commands that are not project-specific.

    Returns the url the user is either logged into or in the case of a user using a personal access token (PAT),
    the default API URL, since PATs don't require login to be used.
    Current examples: auth whoami and some task-history commands that are not project-specific.
    """
    if cli_settings.DIVBASE_API_PAT:
        return cli_settings.DIVBASE_API_URL

    config = load_user_config()
    if not config.logged_in_url:
        raise AuthenticationError("You are not logged in. Please log in with 'divbase-cli auth login [EMAIL]'.")
    return config.logged_in_url


def resolve_project(project_name: str | None) -> ProjectConfig:
    """
    Helper function to resolve the project to use for a CLI command.
    Falls back to the default project set in the user config if not explicitly provided.

    Once the project is resolved a ProjectConfig object is returned,
    which contains the name and API URL of the project.
    """
    config = load_user_config()
    if not project_name:
        project_name = config.default_project
    if not project_name:
        raise ProjectNameNotSpecifiedError()
    return config.project_info(project_name)


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
        return Path(download_dir)
    return Path.cwd()
