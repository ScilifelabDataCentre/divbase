"""
Helper functions for setting up the API for testing.

Creates some users and projects in the test database instance.

TODO: When CLI is more mature, use CLI commands to do this setup instead of direct API calls.
"""

import httpx

BASE_URL = "http://localhost:8001/api"

ADMIN_CREDENTIALS = {"email": "admin@divbase.com", "password": "badpassword"}

TEST_USERS = {
    "read user": {"email": "read@divbase.se", "password": "badpassword"},
    "edit user": {"email": "edit@divbase.se", "password": "badpassword"},
    "manage user": {"email": "manage@divbase.se", "password": "badpassword"},
}

# TODO - bucket and project names should not match to be closer to reality.
TEST_PROJECTS = [
    {
        "name": "project1",
        "description": "First test project",
        "bucket_name": "project1",
        "storage_quota_bytes": 10737418240,
    },
    {
        "name": "project2",
        "description": "Second test project",
        "bucket_name": "project2",
        "storage_quota_bytes": 10737418240,
    },
    {
        "name": "query-project",
        "description": "Third test project",
        "bucket_name": "query-project",
        "storage_quota_bytes": 10737418240,
    },
    {
        "name": "split-scaffold-project",
        "description": "Fourth test project",
        "bucket_name": "split-scaffold-project",
        "storage_quota_bytes": 10737418240,
    },
    {
        "name": "cleaned-project",
        "description": "Fifth test project",
        "bucket_name": "cleaned-project",
        "storage_quota_bytes": 10737418240,
    },
    {
        "name": "empty-project",
        "description": "Sixth test project",
        "bucket_name": "empty-project",
        "storage_quota_bytes": 10737418240,
    },
    {
        "name": "mixed-concat-merge-project",
        "description": "7th test project",
        "bucket_name": "mixed-concat-merge-project",
        "storage_quota_bytes": 10737418240,
    },
]

USER_ROLES = {
    "read@divbase.se": "read",
    "edit@divbase.se": "edit",
    "manage@divbase.se": "manage",
}

ROLE_ASSIGNMENTS = {project["name"]: [(email, role) for email, role in USER_ROLES.items()] for project in TEST_PROJECTS}


def get_admin_access_token() -> str:
    """Get admin access token for authenticated requests."""
    response = httpx.post(
        f"{BASE_URL}/v1/auth/login",
        data={
            "grant_type": "password",
            "username": ADMIN_CREDENTIALS["email"],
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
        )

        user = response.json()
        user_map[user["email"]] = user["id"]

    return user_map


def create_projects(token: str) -> dict[str, int]:
    project_map = {}
    for project_data in TEST_PROJECTS:
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


def setup_api_data() -> None:
    """Create test users and projects in the API."""
    admin_token = get_admin_access_token()
    user_map = create_users(admin_token)
    project_map = create_projects(admin_token)
    assign_project_roles(admin_token, user_map, project_map)
    print("Test users and projects created successfully.")
