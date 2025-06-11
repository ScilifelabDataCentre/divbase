"""
Handles the user's configuration file for the divbase_tools package.
Not to be confused with the bucket versioning configuration.

User configuration is stored in a local file.
By default the config will be stored at: "~/.config/.divbase_tools.yaml"
"""

from pathlib import Path

import yaml

DEFAULT_CONFIG_PATH = Path.home() / ".config" / ".divbase_tools.yaml"


TEMPLATE_USER_CONFIG = {
    "buckets": [],
    "default_bucket": None,
    "DivBase_Access_Key_Env_Name": "DIVBASE_ACCESS_KEY",
    "DivBase_Secret_Key_Env_Name": "DIVBASE_SECRET_KEY",
}


def create_user_config(config_path: Path = DEFAULT_CONFIG_PATH) -> None:
    """Create a user configuration file at the specified path."""
    if config_path.exists():
        raise FileExistsError(f"Config file already exists at {config_path}.")
    config_path.parent.mkdir(parents=False, exist_ok=True)
    with open(config_path, "w") as file:
        yaml.safe_dump(TEMPLATE_USER_CONFIG, file, sort_keys=False)


def add_bucket_to_config(bucket_name: str, config_path: Path, is_default: bool) -> None:
    """
    Add a new bucket to the user configuration file.
    If the file does not exist, it will be created.
    """
    config = load_user_config(config_path)

    if bucket_name not in config["buckets"]:
        config["buckets"].append(bucket_name)

    if is_default:
        config["default_bucket"] = bucket_name

    save_user_config(config, config_path)


def get_default_bucket(config_path: Path) -> str | None:
    """
    Get the default bucket from the user configuration file.
    If no default bucket is set, return None.
    """
    config = load_user_config(config_path)
    return config.get("default_bucket", None)


def set_default_bucket(bucket_name: str, config_path: Path) -> None:
    """
    Set the default bucket in the user configuration file.
    If the file does not exist, it will be created.
    """
    config = load_user_config(config_path)
    config["default_bucket"] = bucket_name
    save_user_config(config, config_path)


def load_user_config(config_path: Path) -> dict:
    with open(config_path, "r") as file:
        return yaml.safe_load(file)


def save_user_config(config: dict, config_path: Path) -> None:
    config_path.parent.mkdir(parents=False, exist_ok=True)
    with open(config_path, "w") as file:
        yaml.safe_dump(config, file, sort_keys=False)
