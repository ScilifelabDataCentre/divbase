"""
A helper script to set up some buckets and add them to user config for DivBase local development.

### This script will:
1. Create some buckets in Minio and populate them with test data from the fixtures directory.
2. Create a user config file for you and add these buckets to it.

### Requirements

1. Have the local docker compose stack running and healthy:

docker compose -f docker/divbase_compose.yaml up -d

# to check if all healthy
docker compose -f docker/divbase_compose.yaml ps

# You could consider this one liner to first clean (DELETE ALL VOLUMES) and start fresh:
docker compose -f docker/divbase_compose.yaml down -v && docker compose -f docker/divbase_compose.yaml up --build -d

### Usage:
Once all containers are healthy, you can run this script with:
uv run scripts/local_dev_setup.py # (or python scripts/local_dev_setup.py)
"""

import contextlib
import os
import shlex
import subprocess
from pathlib import Path

import boto3
import httpx

from divbase_cli.cli_exceptions import DivBaseAPIError

MINIO_URL = "http://localhost:9000"
MINIO_FAKE_ACCESS_KEY = "minioadmin"
MINIO_FAKE_SECRET_KEY = "badpassword"

FIXTURES_DIR = Path(__file__).parent.parent / "tests" / "fixtures"

LOCAL_ENV = os.environ.copy()
LOCAL_ENV["DIVBASE_ENV"] = "local"
BASE_URL = "http://localhost:8000/api"

ADMIN_CREDENTIALS = {"username": "admin@divbase.com", "password": "badpassword"}

USERS_TO_CREATE = [
    {"name": "Alice", "email": "alice@example.com", "password": "badpassword"},
    {"name": "Bob", "email": "bob@example.com", "password": "badpassword"},
    {"name": "Charlie", "email": "charlie@example.com", "password": "badpassword"},
    {"name": "Diana", "email": "diana@example.com", "password": "badpassword"},
]

PROJECTS = [
    {
        "name": "local-project-1",
        "description": "First test project for local development",
        "bucket_name": "divbase-local-1",
        "storage_quota_bytes": 10737418240,
        "files": [
            "HOM_20ind_17SNPs_first_10_samples.vcf.gz",
            "HOM_20ind_17SNPs_last_10_samples.vcf.gz",
            "sample_metadata.tsv",
        ],
    },
    {
        "name": "local-project-2",
        "description": "Second test project for local development",
        "bucket_name": "divbase-local-2",
        "storage_quota_bytes": 10737418240,
        "files": ["sample_metadata.tsv"],
    },
    {
        "name": "local-project-3",
        "description": "Third test project for local development",
        "bucket_name": "divbase-local-3",
        "storage_quota_bytes": 10737418240,
        "files": [
            "HOM_20ind_17SNPs.1.vcf.gz",
            "HOM_20ind_17SNPs.13.vcf.gz",
            "HOM_20ind_17SNPs.18.vcf.gz",
            "HOM_20ind_17SNPs.20.vcf.gz",
            "HOM_20ind_17SNPs.21.vcf.gz",
            "HOM_20ind_17SNPs.22.vcf.gz",
            "HOM_20ind_17SNPs.24.vcf.gz",
            "HOM_20ind_17SNPs.4.vcf.gz",
            "HOM_20ind_17SNPs.5.vcf.gz",
            "HOM_20ind_17SNPs.6.vcf.gz",
            "HOM_20ind_17SNPs.7.vcf.gz",
            "HOM_20ind_17SNPs.8.vcf.gz",
            "sample_metadata_HOM_chr_split_version.tsv",
        ],
    },
    {
        "name": "local-project-4",
        "description": "Fourth test project for local development",
        "bucket_name": "divbase-local-4",
        "storage_quota_bytes": 10737418240,
        "files": [
            "HOM_20ind_17SNPs.1.vcf.gz",
            "HOM_20ind_17SNPs.4.vcf.gz",
            "HOM_20ind_17SNPs.21.vcf.gz",
            "HOM_20ind_17SNPs.8_edit_new_sample_names.vcf.gz",
            "HOM_20ind_17SNPs.13_edit_new_sample_names.vcf.gz",
            "HOM_20ind_17SNPs.18_edit_new_sample_names.vcf.gz",
            "HOM_20ind_17SNPs_changed_sample_names.vcf.gz",
            "sample_metadata_HOM_files_that_need_mixed_bcftools_concat_and_merge.tsv",
        ],
    },
]


ROLE_ASSIGNMENTS = {
    "local-project-1": [
        ("alice@example.com", "manage"),
        ("bob@example.com", "edit"),
        ("charlie@example.com", "read"),
    ],
    "local-project-2": [
        ("alice@example.com", "edit"),
        ("bob@example.com", "manage"),
        ("diana@example.com", "read"),
    ],
    "local-project-3": [
        ("charlie@example.com", "manage"),
        ("diana@example.com", "edit"),
        ("alice@example.com", "read"),
    ],
    "local-project-4": [
        ("diana@example.com", "manage"),
        ("charlie@example.com", "edit"),
        ("bob@example.com", "read"),
        ("alice@example.com", "edit"),
    ],
}


def setup_minio_buckets() -> None:
    """Create buckets and enable bucket versioning"""
    s3_client = boto3.client(
        "s3",
        endpoint_url=MINIO_URL,
        aws_access_key_id=MINIO_FAKE_ACCESS_KEY,
        aws_secret_access_key=MINIO_FAKE_SECRET_KEY,
    )

    for project in PROJECTS:
        with contextlib.suppress(s3_client.exceptions.BucketAlreadyOwnedByYou):
            s3_client.create_bucket(Bucket=project["bucket_name"])

        s3_client.put_bucket_versioning(
            Bucket=project["bucket_name"],
            VersioningConfiguration={"Status": "Enabled"},
        )


def get_admin_access_token() -> str:
    """Get admin access token for authenticated requests."""
    response = httpx.post(
        f"{BASE_URL}/v1/auth/login",
        data={
            "grant_type": "password",
            "username": ADMIN_CREDENTIALS["username"],
            "password": ADMIN_CREDENTIALS["password"],
            "scope": "",
            "client_id": "string",
            "client_secret": "",
        },
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        timeout=10.0,
    )
    response.raise_for_status()

    data = response.json()
    return data["access_token"]


def make_authenticated_request(method: str, url: str, token: str, **kwargs) -> httpx.Response:
    """Appends token etc... to make an authenticated request to API."""
    headers = kwargs.get("headers", {})
    headers["Authorization"] = f"Bearer {token}"
    kwargs["headers"] = headers
    kwargs["timeout"] = 10.0

    response = httpx.request(method, url, **kwargs)

    try:
        response.raise_for_status()
    except httpx.HTTPStatusError:
        error_details = response.json().get("detail", "No error details provided")
        error_type = response.json().get("type", "unknown")
        raise DivBaseAPIError(
            error_details=error_details,
            status_code=response.status_code,
            http_method=method,
            url=url,
            error_type=error_type,
        ) from None

    return response


def create_users(token: str) -> dict[str, int]:
    print("Creating test users...")

    user_map = {}
    for user_data in USERS_TO_CREATE:
        response = make_authenticated_request(
            "POST",
            f"{BASE_URL}/v1/admin/users/",
            token,
            json={"name": user_data["name"], "email": user_data["email"], "password": user_data["password"]},
            params={"email_verified": True},
        )

        user = response.json()
        user_map[user["email"]] = user["id"]
        print(f"Created user: {user['name']} ({user['email']})")

    return user_map


def create_projects(token: str) -> dict[str, int]:
    print("Creating test projects...")
    project_map = {}
    for project in PROJECTS:
        project_data = {
            "name": project["name"],
            "description": project["description"],
            "bucket_name": project["bucket_name"],
            "storage_quota_bytes": project["storage_quota_bytes"],
        }
        response = make_authenticated_request("POST", f"{BASE_URL}/v1/admin/projects", token, json=project_data)

        project = response.json()
        project_map[project["name"]] = project["id"]
        print(f"Created project: {project['name']}")

    return project_map


def assign_project_roles(token: str, user_map: dict[str, int], project_map: dict[str, int]) -> None:
    """Assign users to projects with specified roles."""
    print("Assigning users to projects with roles...")

    for project_name, assignments in ROLE_ASSIGNMENTS.items():
        project_id = project_map[project_name]
        print(f"Assigning roles for project: {project_name}")

        for user_email, role in assignments:
            user_id = user_map[user_email]

            make_authenticated_request(
                "POST", f"{BASE_URL}/v1/admin/projects/{project_id}/members/{user_id}", token, params={"role": role}
            )
            print(f"Assigned {user_email} as {role} to {project_name}")

        # assign admin as manager to all projects
        # Hardcoded user_id=1 as FIRST_ADMIN_USER should have that by being created in local_dev_setup.py
        make_authenticated_request(
            "POST", f"{BASE_URL}/v1/admin/projects/{project_id}/members/1", token, params={"role": "manage"}
        )
        print(f"Assigned first admin user as manager to {project_name}")


def create_local_config():
    """Create a local user DivBase config file and add projects to it."""
    command = shlex.split("divbase-cli config create")
    result = subprocess.run(command, check=False, env=LOCAL_ENV, stderr=subprocess.PIPE)
    if result.returncode != 0 and "FileExistsError" not in result.stderr.decode():
        raise subprocess.CalledProcessError(result.returncode, command, output=result.stdout, stderr=result.stderr)

    for project in PROJECTS:
        command = shlex.split(f"divbase-cli config add {project['name']}")
        subprocess.run(command, check=True, env=LOCAL_ENV)

    default_project = PROJECTS[0]["name"]
    command = shlex.split(f"divbase-cli config set-default {default_project}")
    subprocess.run(command, check=True, env=LOCAL_ENV)


def login_to_divbase():
    """Log in to DivBase using the local environment."""
    command = shlex.split("divbase-cli auth login admin@divbase.com --password badpassword")
    subprocess.run(command, check=True, env=LOCAL_ENV)


def upload_files_to_buckets():
    """Upload files to each created project's storage bucket."""
    for project in PROJECTS:
        project_name = project["name"]
        files = project["files"]
        file_list = [str(FIXTURES_DIR / file) for file in files]
        files_to_upload = " ".join(file_list)
        command = shlex.split(f"divbase-cli files upload --project {project_name} {files_to_upload}")
        subprocess.run(command, check=True, env=LOCAL_ENV)


def create_first_project_version():
    """Add a user defined project version to each project after initial file upload."""
    for project in PROJECTS:
        command = shlex.split(
            f"divbase-cli version add v0.1.0 --description 'add initial data sets' --project {project['name']}"
        )
        subprocess.run(command, check=True, env=LOCAL_ENV)


if __name__ == "__main__":
    setup_minio_buckets()

    admin_token = get_admin_access_token()
    user_map = create_users(admin_token)
    project_map = create_projects(admin_token)
    assign_project_roles(admin_token, user_map, project_map)

    create_local_config()
    login_to_divbase()
    upload_files_to_buckets()
    create_first_project_version()

    print("Setup completed successfully!")
    print("\nTest Users:")
    for user_data in USERS_TO_CREATE:
        print(f"   {user_data['email']} (password: {user_data['password']})")

    print("\nTest Projects:")
    for project_name in project_map:
        print(f"   {project_name}")
