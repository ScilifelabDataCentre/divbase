"""
Functions that resolve for the CLI commands things like:
    - which project (and associated storage bucket) to use
    - which download directory to use
    - which DivBase API/S3 URL to use
Based on provider user input and their config file.
"""

from pathlib import Path

from divbase_tools.exceptions import ProjectNameNotSpecifiedError
from divbase_tools.user_config import ProjectConfig, load_user_config


def resolve_project(project_name: str | None, config_path: Path) -> ProjectConfig:
    """
    Helper function to resolve the project to use for a CLI command.
    Falls back to the default project set in the user config if not explicitly provided.

    Once the project is resolved a ProjectConfig object is returned,
    which contains name and URLs (S3+API) for the project.
    """
    config = load_user_config(config_path)
    if not project_name:
        project_name = config.default_project
    if not project_name:
        raise ProjectNameNotSpecifiedError(config_path=config_path)
    return config.project_info(project_name)


def resolve_divbase_api_url(url: str | None, config_path: Path) -> str:
    """
    Helper function to resolve the DivBase API URL to use for a CLI command.

    If not provided by the user, it will take the default project's API URL from the user config.
    Otherwise raise an error.
    """
    if url:
        return url

    config = load_user_config(config_path=config_path)

    if not config.default_project:
        raise ValueError(
            "No default project is set in your user config. "
            "Please set a default project or specify the API URL using the --url option."
        )

    project_config = config.project_info(name=config.default_project)
    return project_config.divbase_url


def resolve_download_dir(download_dir: str | None, config_path: Path) -> Path:
    """
    Helper function to resolve the download directory to use for a CLI command involving downloading files.

    Priority given to `download_dir` argument, then if a default is set in the user config.
    Note: "." or None should default to the current working directory.
    """
    if not download_dir:
        config = load_user_config(config_path)
        download_dir = config.default_download_dir

    if download_dir and download_dir != ".":
        return Path(download_dir)
    return Path.cwd()
