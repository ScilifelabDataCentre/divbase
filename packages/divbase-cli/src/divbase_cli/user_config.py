"""
Handles the user's configuration file for the divbase-cli package.
User configuration is stored in a local file.
By default the config will be stored at: "~/.config/.divbase_tools.yaml" # TODO - change name?

Not to be confused with the bucket versioning files, which is stored in the project's bucket
and versions the state of all files in the bucket at given timestamps.
"""

import warnings
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from divbase_cli.cli_config import cli_settings
from divbase_lib.exceptions import ProjectNotInConfigError


@dataclass
class ProjectConfig:
    """
    Config of a single project.

    TODO - consider adding validation of URLs and names based on known constraints.
    """

    name: str
    divbase_url: str
    s3_url: str

    @property
    def bucket_name(self) -> str:
        """
        The name of the storage project associated with this project.

        Note: As of now it is the same, but in the future it might be different, but hopefully can be derived from the project name directly.
        """
        return self.name


@dataclass
class UserConfig:
    """Overall user configuration"""

    config_path: Path
    logged_in_url: str | None = None  # URL of the divbase server the user is currently logged in to, if any.
    projects: list[ProjectConfig] = field(default_factory=list)
    default_project: str | None = None
    # default_download_dir is a string (rather than Path) for easier loading/saving.
    default_download_dir: str | None = None

    @property
    def all_project_names(self) -> list[str]:
        """List of all project names in the user's config"""
        return [project.name for project in self.projects]

    def dump_config(self) -> None:
        """
        Dump the user configuration to the specified output path
        Don't include the config_path in the dumped file.
        """
        config_dict = {
            "logged_in_url": self.logged_in_url,
            "projects": [project.__dict__ for project in self.projects],
            "default_project": self.default_project,
            "default_download_dir": self.default_download_dir,
        }

        with open(self.config_path, "w") as file:
            yaml.safe_dump(config_dict, file, sort_keys=False)

    def add_project(self, name: str, divbase_url: str, s3_url: str, is_default: bool) -> ProjectConfig:
        """
        Add a new project to the user configuration file.
        If the configuration file does not exist, it will be created.
        """
        new_project = ProjectConfig(name=name, divbase_url=divbase_url, s3_url=s3_url)

        if new_project.name in self.all_project_names:
            warnings.warn(
                f"""The project: '{new_project.name}' already existed in your configuration file. It will be replaced with your new params""",
                stacklevel=2,
            )

            self.projects = [project for project in self.projects if project.name != new_project.name]

        self.projects.append(new_project)

        if is_default:
            self.default_project = new_project.name

        self.dump_config()
        return new_project

    def set_default_project(self, name: str) -> str:
        """
        Set the default project in the user configuration file.
        """
        if name not in self.all_project_names:
            raise ValueError(
                f"The project name specified: '{name}', was not found in your config file, please add it first."
            )
        self.default_project = name
        self.dump_config()
        return self.default_project

    def remove_project(self, name: str) -> str | None:
        """
        Remove a project from the user configuration file.
        Returns the project name if it was removed successfully, otherwise None.
        """
        if name not in self.all_project_names:
            return

        self.projects = [project for project in self.projects if project.name != name]

        if self.default_project == name:
            self.default_project = None

        self.dump_config()

        return name

    def set_default_download_dir(self, download_dir: str) -> str:
        """
        Set the default download directory in the user configuration file.
        """
        self.default_download_dir = download_dir
        self.dump_config()
        return self.default_download_dir

    def project_info(self, name: str) -> ProjectConfig:
        """
        Get the project configuration given the project name.
        """
        for project in self.projects:
            if project.name == name:
                return project
        raise ProjectNotInConfigError(config_path=self.config_path, project_name=name)

    def set_logged_in_url(self, url: str | None) -> None:
        """
        Set the URL of the DivBase server the user is currently logged in to.
        Used during login/logout to track the state.
        """
        self.logged_in_url = url
        self.dump_config()


def load_user_config(config_path: Path = cli_settings.CONFIG_PATH) -> UserConfig:
    """Helper function to load the user config file"""
    with open(config_path, "r") as file:
        config_contents = yaml.safe_load(file)

    projects = [ProjectConfig(**project) for project in config_contents.get("projects", [])]

    return UserConfig(
        logged_in_url=config_contents.get("logged_in_url"),
        config_path=config_path,
        projects=projects,
        default_project=config_contents.get("default_project"),
        default_download_dir=config_contents.get("default_download_dir"),
    )


def create_user_config(config_path: Path) -> None:
    """Create a user configuration file at the specified path."""
    if config_path.exists():
        raise FileExistsError(f"Config file already exists at {config_path}.")
    config_path.parent.mkdir(parents=False, exist_ok=True)

    user_config = UserConfig(config_path=config_path, projects=[])
    user_config.dump_config()
