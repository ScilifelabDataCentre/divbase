import contextlib
import os
import shlex
import subprocess
from pathlib import Path

import boto3
import requests
import yaml

MINIO_URL = "http://localhost:9000"
MINIO_FAKE_ACCESS_KEY = "minioadmin"
MINIO_FAKE_SECRET_KEY = "badpassword"

FIXTURES_DIR = Path(__file__).parent.parent / "tests" / "fixtures"

LOCAL_ENV = os.environ.copy()
LOCAL_ENV["DIVBASE_ENV"] = "local"

DIVBASE_API_URL = "http://localhost:8001/api/v1"
DIVBASE_ADMIN_EMAIL = "admin@divbase.com"
DIVBASE_ADMIN_PASSWORD = "badpassword"


def get_divbase_access_token() -> str:
    """Get an access token for the DivBase API."""
    response = requests.post(
        f"{DIVBASE_API_URL}/../auth/login",
        data={"username": DIVBASE_ADMIN_EMAIL, "password": DIVBASE_ADMIN_PASSWORD},
    )
    response.raise_for_status()
    return response.json()["access_token"]


def ensure_project_exists(project_name: str):
    """Create a Minio bucket and DivBase project config entry if they do not already exist."""

    print(f"\nChecking if there is a bucket named {project_name} in the local MinIO container...")

    s3_client = boto3.client(
        "s3",
        endpoint_url=MINIO_URL,
        aws_access_key_id=MINIO_FAKE_ACCESS_KEY,
        aws_secret_access_key=MINIO_FAKE_SECRET_KEY,
    )

    existing_buckets = [b["Name"] for b in s3_client.list_buckets().get("Buckets", [])]
    if project_name in existing_buckets:
        print(f"Bucket '{project_name}' already exists.")
    else:
        print(f"Bucket '{project_name}' does not exist. Creating it...")
        with contextlib.suppress(s3_client.exceptions.BucketAlreadyOwnedByYou):
            s3_client.create_bucket(Bucket=project_name)
        print(f"Bucket '{project_name}' created.")

    s3_client.put_bucket_versioning(
        Bucket=project_name,
        VersioningConfiguration={"Status": "Enabled"},
    )


def ensure_project_in_config(project_name: str) -> None:
    """Ensure that the project is listed in the DivBase config."""

    config_path = Path.home() / ".config" / ".divbase_tools.yaml"

    print(f"\nChecking if {project_name} exists in config...")

    if project_in_config(project_name, config_path):
        print(f"Project '{project_name}' already exists in config.")
    else:
        print(f"Project '{project_name}' does not exist in config. Adding it...")
        command = shlex.split(f"divbase-cli config add {project_name}")
        subprocess.run(command, check=True, env=LOCAL_ENV)


def project_in_config(project_name: str, config_path: Path) -> bool:
    if not config_path.exists():
        return False
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    projects = config.get("projects", [])
    return any(p.get("name") == project_name for p in projects)
