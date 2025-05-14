"""
Manage versions of the entire bucket

This script is interacted with via argparse. 

Versioning specified by updating a yaml file in the bucket called ".bucket_versions.yaml".

Versioning done by setting the version according to the current timestamp.  
"""
from dataclasses import dataclass, field
import yaml
import argparse
from pathlib import Path
from datetime import datetime, timezone

BUCKET_VERSION_FILE = ".bucket_versions.yaml"

@dataclass
class BucketVersionManager:
    base_dir: Path
    # these 2 are not defined till post_init method
    version_file: Path = field(init=False)
    version_info: dict = field(init=False)

    def __post_init__(self):
        self.base_dir = Path(self.base_dir).resolve()
        self.version_file = self.base_dir / BUCKET_VERSION_FILE
        self.version_info = {}

        if not self.version_file.exists():
            return

        with self.version_file.open("r") as f:
            self.version_info = yaml.safe_load(f)
            print(f"Loaded version info from file: {self.version_file.resolve()}")

    def create_metadata_file(self):
        """Create the initial metadata file with a default version."""
        if self.version_file.exists():
            print(f"Version file already exists at {self.version_file.resolve()}.")
            return

        timestamp = self._create_timestamp()
        content = {
            "versions": {
                "v0.1.0": {
                    "timestamp": timestamp,
                    "description": "Initial version"
                }
            }
        }

        self._write_to_file(content)
        print(f"Metadata file created at {self.version_file.resolve()}.")

    def add_version(self, name: str, description: str):
        """Add a new version to the metadata file."""
        if not self.version_file.exists():
            print("No version file found in this directory, ensure you're in the correct place.")
            return
        
        version_data = self.version_info

        # Add the new version
        timestamp = self._create_timestamp()
        version_data["versions"][name] = {
            "timestamp": timestamp,
            "description": description
        }

        self._write_to_file(version_data)
        # only update the version_info if the file was successfully written
        self.version_info = version_data
        print(f"Version {name} added to version file at {self.version_file.resolve()}.")

    def list_versions(self):

        # TODO - think about this logic of 
        # checking for both file existence and version info in the class
        if not self.version_file.exists():
            print("No version file found in this directory, ensure you're in the correct place.")
            return
        if not self.version_info:
            print("No versioning information present as of yet")
            return

        print("Versions in the metadata file:")
        for version, details in self.version_info["versions"].items():
            print(f"- {version}: {details['timestamp']} ({details['description']})")

    def _create_timestamp(self) -> str:
        return datetime.now(tz=timezone.utc).isoformat()

    def _read_from_file(self) -> dict:
        with self.version_file.open("r") as f:
            return yaml.safe_load(f)

    def _write_to_file(self, content: dict):
        with self.version_file.open("w") as f:
            yaml.safe_dump(content, f, default_flow_style=False)


def create_command(args):
    manager = BucketVersionManager(args.base_dir)
    manager.create_metadata_file()

def add_command(args):
    manager = BucketVersionManager(args.base_dir)
    manager.add_version(args.name, args.description)

def list_command(args):
    manager = BucketVersionManager(args.base_dir)
    manager.list_versions()

def delete_command(args):
    pass



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Manage bucket versions")

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("--base_dir", default=".", help="Base directory for the version metadata file.")

    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: create
    create_parser = subparsers.add_parser("create", parents=[parent_parser], help="Create the versioning metadata file.")
    create_parser.set_defaults(func=create_command)

    # Subcommand: add
    add_parser = subparsers.add_parser("add", parents=[parent_parser], help="Add a new version.")
    add_parser.add_argument("--name", required=True, help="Name of the version (e.g., semantic version).")
    add_parser.add_argument("--description", required=False, help="Optional description of the version.")
    add_parser.set_defaults(func=add_command)

    # Subcommand: list
    list_parser = subparsers.add_parser("list", parents=[parent_parser], help="List all versions.")
    list_parser.set_defaults(func=list_command)

    args = parser.parse_args()
    args.func(args)
