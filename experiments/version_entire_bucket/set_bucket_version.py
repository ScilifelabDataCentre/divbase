"""
Manage versions of the entire bucket

This script is interacted with via argparse. 

Versioning specified by updating a yaml file in the bucket called ".bucket_versions.yaml".

Versioning done by setting the version according to the current timestamp.  
"""
import yaml
import argparse
from pathlib import Path
from datetime import datetime, timezone

BUCKET_VERSION_FILE = ".bucket_versions.yaml"


def create_metadata_file(base_dir: str):
    version_file = Path(base_dir) / BUCKET_VERSION_FILE
    if version_file.exists():
        print(f"Version file already exists at {version_file.resolve()}.")
        return
    
    timestamp = create_timestamp()

    content = {
        "versions": {
           "v0.1.0" : {
                "timestamp": timestamp,
                "description": "Initial version"
            }
        }
    }

    with version_file.open("w") as f:
        yaml.safe_dump(content, f)
    print(f"Metadata file created at {version_file.resolve()}.")


def add_version(base_dir: str, name: str, description: str):
    """
    Add a new version to the metadata file.
    Question - should this be auto added to the bucket? 

    # TODO - add check name does not exist

    """
    version_file = Path(base_dir) / BUCKET_VERSION_FILE
    if not version_file.exists():
        print("No version file found in this directory, ensure you're in the correct place")
        return
    print(f"Adding version {name} with description: {description}")

    with version_file.open("r") as f:
        version_data = yaml.safe_load(f)

    timestamp = create_timestamp()
        
    version_data["versions"][name] = {
        "timestamp": timestamp,
        "description": description
    }

    
    with version_file.open("w") as f:
        yaml.safe_dump(version_data, f)
    print(f"Version {name} added to version file at {version_file.resolve()}.")



def list_versions(base_dir: str):
    print("Listing all versions...")



def validate_version_file():
    pass


def create_timestamp() -> str:
    return datetime.now(tz=timezone.utc).isoformat()





if __name__ == "__main__":

    create_timestamp()

    parser = argparse.ArgumentParser(description="Manage bucket versions")

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("--base_dir", default=".", help="Base directory for the version metadata file.")

    subparsers = parser.add_subparsers(dest="command", required=True)

    create_parser = subparsers.add_parser("create", parents=[parent_parser], help="Create the versioning metadata file.")
    create_parser.set_defaults(func=lambda args: create_metadata_file(args.base_dir))

    add_parser = subparsers.add_parser("add", parents=[parent_parser], help="Add a new version.")
    add_parser.add_argument("--name", required=True, help="Name of the version (e.g., semantic version).")
    add_parser.add_argument("--description", required=False, help="Optional description of the version.")
    add_parser.set_defaults(func=lambda args: add_version(args.base_dir, args.name, args.description))

    list_parser = subparsers.add_parser("list", parents=[parent_parser], help="List all versions.")
    list_parser.set_defaults(func=lambda args: print("List all versions."))


    args = parser.parse_args()
    args.func(args)

