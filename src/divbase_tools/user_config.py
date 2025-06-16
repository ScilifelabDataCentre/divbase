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
    """
    config = load_user_config(config_path)

    if bucket_name not in config["buckets"]:
        config["buckets"].append(bucket_name)
    config["default_bucket"] = bucket_name

    save_user_config(config, config_path)


def remove_bucket_from_config(bucket_name: str, config_path: Path) -> str:
    """
    Remove a bucket from the user configuration file.
    Returns the bucket name if it was removed successfully.
    """
    config = load_user_config(config_path)

    if bucket_name in config["buckets"]:
        config["buckets"].remove(bucket_name)
    else:
        raise ValueError(f"The bucket specified: '{bucket_name}', was not found in your config file.")

    if config.get("default_bucket") == bucket_name:
        config["default_bucket"] = None

    save_user_config(config, config_path)


def load_user_config(config_path: Path) -> dict:
    with open(config_path, "r") as file:
        return yaml.safe_load(file)


def save_user_config(config: dict, config_path: Path) -> None:
    config_path.parent.mkdir(parents=False, exist_ok=True)
    with open(config_path, "w") as file:
        yaml.safe_dump(config, file, sort_keys=False)
