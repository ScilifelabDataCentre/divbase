"""
Handles the user's configuration file for the divbase_tools package.
Not to be confused with the bucket versioning configuration.

User configuration is stored in a local file.
By default the config will be stored at: "~/.config/.divbase_tools.yaml"
"""

import warnings
from dataclasses import dataclass, field
from pathlib import Path

import yaml


@dataclass
class BucketConfig:
    """
    Config of a single bucket.

    TODO - consider adding validation of URLs and names based on known constaints.
    """

    name: str
    divbase_url: str
    s3_url: str


@dataclass
class UserConfig:
    """
    Overall user configuration

    """

    config_path: Path
    buckets: list[BucketConfig] = field(default_factory=list)
    default_bucket: str | None = None
    # default_download_dir is a string (rather than Path) for easier loading/saving.
    default_download_dir: str | None = None

    @property
    def all_bucket_names(self) -> list[str]:
        """List of all bucket names in the user's config"""
        return [bucket.name for bucket in self.buckets]

    def dump_config(self) -> None:
        """
        Dump the user configuration to the specified output path
        Don't include the config_path in the dumped file.
        """
        config_dict = {
            "buckets": [bucket.__dict__ for bucket in self.buckets],
            "default_bucket": self.default_bucket,
            "default_download_dir": self.default_download_dir,
        }

        with open(self.config_path, "w") as file:
            yaml.safe_dump(config_dict, file, sort_keys=False)

    def add_bucket(self, name: str, divbase_url: str, s3_url: str, is_default: bool) -> BucketConfig:
        """
        Add a new bucket to the user configuration file.
        If the configuration file does not exist, it will be created.
        """
        new_bucket = BucketConfig(name=name, divbase_url=divbase_url, s3_url=s3_url)

        if new_bucket.name in self.all_bucket_names:
            warnings.warn(
                f"""The bucket: '{new_bucket.name}' already existed in your configuration file. It will be replaced with your new params""",
                stacklevel=2,
            )

            self.buckets = [bucket for bucket in self.buckets if bucket.name != new_bucket.name]

        self.buckets.append(new_bucket)

        if is_default:
            self.default_bucket = new_bucket.name

        self.dump_config()
        return new_bucket

    def set_default_bucket(self, name: str) -> str:
        """
        Set the default bucket in the user configuration file.
        """
        if name not in self.all_bucket_names:
            raise ValueError(
                f"The bucket name specified: '{name}', was not found in your config file, please add it first."
            )
        self.default_bucket = name
        self.dump_config()
        return self.default_bucket

    def remove_bucket(self, name: str) -> str | None:
        """
        Remove a bucket from the user configuration file.
        Returns the bucket name if it was removed successfully, otherwise None.
        """
        if name not in self.all_bucket_names:
            return

        self.buckets = [bucket for bucket in self.buckets if bucket.name != name]

        if self.default_bucket == name:
            self.default_bucket = None

        self.dump_config()

        return name

    def set_default_download_dir(self, download_dir: str) -> str:
        """
        Set the default download directory in the user configuration file.
        """
        self.default_download_dir = download_dir
        self.dump_config()
        return self.default_download_dir

    def bucket_info(self, name: str) -> BucketConfig | None:
        """
        Get the bucket configuration given the bucket name.
        """
        for bucket in self.buckets:
            if bucket.name == name:
                return bucket


def load_user_config(config_path: Path) -> UserConfig:
    """Helper function to load the user config file"""
    with open(config_path, "r") as file:
        config_contents = yaml.safe_load(file)

    buckets = [BucketConfig(**bucket) for bucket in config_contents.get("buckets", [])]

    return UserConfig(
        config_path=config_path,
        buckets=buckets,
        default_bucket=config_contents.get("default_bucket"),
        default_download_dir=config_contents.get("default_download_dir"),
    )


def create_user_config(config_path: Path) -> None:
    """Create a user configuration file at the specified path."""
    if config_path.exists():
        raise FileExistsError(f"Config file already exists at {config_path}.")
    config_path.parent.mkdir(parents=False, exist_ok=True)

    user_config = UserConfig(config_path=config_path, buckets=[])
    user_config.dump_config()
