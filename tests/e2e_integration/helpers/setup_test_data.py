"""
Provides the docker testing stack with test data.

Creates test projects and users with different access roles to those projects
Assigns each project an S3 (MinIO) bucket and populates those buckets with different test files
"""

from pathlib import Path

import boto3
import httpx

BASE_URL = "http://localhost:8001/api"

MINIO_URL = "http://localhost:9002"
MINIO_FAKE_ACCESS_KEY = "minioadmin"  # we need to use admin credentials here to create buckets
MINIO_FAKE_SECRET_KEY = "badpassword"

FIXTURES_DIR = Path(__file__).parent.parent.parent / "fixtures"

API_ADMIN_CREDENTIALS = {"email": "admin@divbase.com", "password": "badpassword"}

TEST_USERS = {
    "read user": {"email": "read@divbase.se", "password": "badpassword"},
    "edit user": {"email": "edit@divbase.se", "password": "badpassword"},
    "manage user": {"email": "manage@divbase.se", "password": "badpassword"},
    "edit user query-project only": {"email": "edit_query_project_only@divbase.se", "password": "badpassword"},
    "manage user query-project only": {"email": "manage_query_project_only@divbase.se", "password": "badpassword"},
}


TEST_PROJECTS = {
    "project1": {
        "description": "First test project",
        "bucket_name": "divbase-local-1",
        "storage_quota_bytes": 10737418240,
        "files": ["file1.txt", "file2.txt"],
    },
    "project2": {
        "description": "Second test project",
        "bucket_name": "divbase-local-2",
        "storage_quota_bytes": 10737418240,
        "files": ["file1.txt"],
    },
    "query-project": {
        "description": "Third test project",
        "bucket_name": "divbase-local-query-project",
        "storage_quota_bytes": 10737418240,
        "files": [
            "HOM_20ind_17SNPs_first_10_samples.vcf.gz",
            "HOM_20ind_17SNPs_last_10_samples.vcf.gz",
            "sample_metadata.tsv",
        ],
    },
    "split-scaffold-project": {
        "description": "Fourth test project",
        "bucket_name": "divbase-local-split-scaffold-project",
        "storage_quota_bytes": 10737418240,
        "files": [
            "HOM_20ind_17SNPs.1.vcf.gz",
            "HOM_20ind_17SNPs.4.vcf.gz",
            "HOM_20ind_17SNPs.5.vcf.gz",
            "HOM_20ind_17SNPs.6.vcf.gz",
            "HOM_20ind_17SNPs.7.vcf.gz",
            "HOM_20ind_17SNPs.8.vcf.gz",
            "HOM_20ind_17SNPs.13.vcf.gz",
            "HOM_20ind_17SNPs.18.vcf.gz",
            "HOM_20ind_17SNPs.20.vcf.gz",
            "HOM_20ind_17SNPs.21.vcf.gz",
            "HOM_20ind_17SNPs.22.vcf.gz",
            "HOM_20ind_17SNPs.24.vcf.gz",
            "sample_metadata_HOM_chr_split_version.tsv",
        ],
    },
    "cleaned-project": {
        "description": "Fifth test project",
        "bucket_name": "divbase-local-cleaned-project",
        "storage_quota_bytes": 10737418240,
        "files": [],  # this project's bucket is always cleaned before a test
    },
    "empty-project": {
        "description": "Sixth test project",
        "bucket_name": "divbase-local-empty-project",
        "storage_quota_bytes": 10737418240,
        "files": [],
    },
    "mixed-concat-merge-project": {
        "description": "7th test project",
        "bucket_name": "divbase-local-mixed-concat-merge-project",
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
}

USER_ROLES = {
    "read@divbase.se": "read",
    "edit@divbase.se": "edit",
    "manage@divbase.se": "manage",
}

ROLE_ASSIGNMENTS = {project: [(email, role) for email, role in USER_ROLES.items()] for project in TEST_PROJECTS}

# Add two special users that only belong to query-project
ROLE_ASSIGNMENTS["query-project"].extend(
    [
        ("edit_query_project_only@divbase.se", "edit"),
        ("manage_query_project_only@divbase.se", "manage"),
    ]
)

_PROJECT_MAP_CACHE: dict[str, int] | None = None


def setup_minio_data() -> None:
    """Create test buckets and add some files to them in Minio."""
    s3_client = boto3.client(
        "s3",
        endpoint_url=MINIO_URL,
        aws_access_key_id=MINIO_FAKE_ACCESS_KEY,
        aws_secret_access_key=MINIO_FAKE_SECRET_KEY,
    )

    for proj_details in TEST_PROJECTS.values():
        bucket_name = proj_details["bucket_name"]
        s3_client.create_bucket(Bucket=bucket_name)
        s3_client.put_bucket_versioning(
            Bucket=bucket_name,
            VersioningConfiguration={"Status": "Enabled"},
        )
        files = proj_details["files"]
        for file in files:
            s3_client.upload_file(Filename=str(FIXTURES_DIR / file), Bucket=bucket_name, Key=file)


def get_admin_access_token() -> str:
    """Get admin access token for authenticated requests."""
    response = httpx.post(
        f"{BASE_URL}/v1/auth/login",
        data={
            "grant_type": "password",
            "username": API_ADMIN_CREDENTIALS["email"],
            "password": API_ADMIN_CREDENTIALS["password"],
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
    response.raise_for_status()
    return response


def create_users(token: str) -> dict[str, int]:
    user_map = {}
    for name, creds in TEST_USERS.items():
        response = make_authenticated_request(
            "POST",
            f"{BASE_URL}/v1/admin/users/",
            token,
            json={"name": name, "email": creds["email"], "password": creds["password"]},
            params={"email_verified": True},
        )

        user = response.json()
        user_map[user["email"]] = user["id"]

    return user_map


def create_projects(token: str) -> dict[str, int]:
    project_map = {}
    for project in TEST_PROJECTS:
        project_data = {
            "name": project,
            "description": TEST_PROJECTS[project]["description"],
            "bucket_name": TEST_PROJECTS[project]["bucket_name"],
            "storage_quota_bytes": TEST_PROJECTS[project]["storage_quota_bytes"],
        }
        response = make_authenticated_request("POST", f"{BASE_URL}/v1/admin/projects", token, json=project_data)

        project = response.json()
        project_map[project["name"]] = project["id"]

    return project_map


def assign_project_roles(token: str, user_map: dict[str, int], project_map: dict[str, int]) -> None:
    """Assign users to projects with specified roles."""
    for project_name, assignments in ROLE_ASSIGNMENTS.items():
        project_id = project_map[project_name]

        for user_email, role in assignments:
            user_id = user_map[user_email]

            make_authenticated_request(
                "POST", f"{BASE_URL}/v1/admin/projects/{project_id}/members/{user_id}", token, params={"role": role}
            )


def get_project_map() -> dict[str, int]:
    """Get the cached project map created during setup."""
    if _PROJECT_MAP_CACHE is None:
        raise RuntimeError("Project map not available. API setup has not been run yet.")
    return _PROJECT_MAP_CACHE


def setup_test_data() -> None:
    """Call all setup functions in the correct order to give the testing stack some test data."""
    setup_minio_data()
    admin_token = get_admin_access_token()
    user_map = create_users(admin_token)
    project_map = create_projects(admin_token)
    assign_project_roles(admin_token, user_map, project_map)
    global _PROJECT_MAP_CACHE
    _PROJECT_MAP_CACHE = project_map
    print("Test users and projects created successfully.")
