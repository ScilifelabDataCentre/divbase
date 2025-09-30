"""
Local dev convenience script.

Separated from local_dev_setup.py until happy with approach.
Will likely convert these to divbase0-cli command first.

Creates test users and projects with various role assignments for development.
Assumes the DivBase stack is already running locally on http://localhost:8000
"""

import httpx

BASE_URL = "http://localhost:8000/api/v1"

ADMIN_CREDENTIALS = {"username": "admin@divbase.com", "password": "badpassword"}

USERS_TO_CREATE = [
    {"name": "Alice", "email": "alice@example.com", "password": "badpassword"},
    {"name": "Bob", "email": "bob@example.com", "password": "badpassword"},
    {"name": "Charlie", "email": "charlie@example.com", "password": "badpassword"},
    {"name": "Diana", "email": "diana@example.com", "password": "badpassword"},
]

PROJECTS_TO_CREATE = [
    {
        "name": "local-squirrel-1",
        "description": "First test project for local development",
        "bucket_name": "local-project-1",
        "storage_quota_bytes": 10737418240,
    },
    {
        "name": "local-mongoose-2",
        "description": "Second test project for local development",
        "bucket_name": "local-project-2",
        "storage_quota_bytes": 10737418240,
    },
    {
        "name": "local-salmon-3",
        "description": "Third test project for local development",
        "bucket_name": "local-project-3",
        "storage_quota_bytes": 10737418240,
    },
    {
        "name": "local-badger-4",
        "description": "Fourth test project for local development",
        "bucket_name": "local-project-4",
        "storage_quota_bytes": 10737418240,
    },
]

ROLE_ASSIGNMENTS = {
    "local-squirrel-1": [
        ("alice@example.com", "manage"),
        ("bob@example.com", "edit"),
        ("charlie@example.com", "read"),
    ],
    "local-mongoose-2": [
        ("alice@example.com", "edit"),
        ("bob@example.com", "manage"),
        ("diana@example.com", "read"),
    ],
    "local-salmon-3": [
        ("charlie@example.com", "manage"),
        ("diana@example.com", "edit"),
        ("alice@example.com", "read"),
    ],
    "local-badger-4": [
        ("diana@example.com", "manage"),
        ("charlie@example.com", "edit"),
        ("bob@example.com", "read"),
        ("alice@example.com", "edit"),
    ],
}


def get_admin_access_token() -> str:
    """Get admin access token for authenticated requests."""
    print("Getting admin access token...")

    response = httpx.post(
        f"{BASE_URL}/auth/login",
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
    print(f"Successfully authenticated as {ADMIN_CREDENTIALS['username']}")
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
    print("Creating test users...")

    user_map = {}
    for user_data in USERS_TO_CREATE:
        response = make_authenticated_request(
            "POST",
            f"{BASE_URL}/admin/users/",
            token,
            json={"name": user_data["name"], "email": user_data["email"], "password": user_data["password"]},
        )

        user = response.json()
        user_map[user["email"]] = user["id"]
        print(f"Created user: {user['name']} ({user['email']})")

    return user_map


def create_projects(token: str) -> dict[str, int]:
    print("Creating test projects...")
    project_map = {}
    for project_data in PROJECTS_TO_CREATE:
        response = make_authenticated_request("POST", f"{BASE_URL}/admin/projects", token, json=project_data)

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
                "POST", f"{BASE_URL}/admin/projects/{project_id}/members/{user_id}", token, params={"role": role}
            )
            print(f"Assigned {user_email} as {role} to {project_name}")

        # assign admin as manager to all projects
        # Hardcoded user_id=1 as FIRST_ADMIN_USER should have that by being created in local_dev_setup.py
        make_authenticated_request(
            "POST", f"{BASE_URL}/admin/projects/{project_id}/members/1", token, params={"role": "manage"}
        )
        print(f"Assigned first admin user as manager to {project_name}")


if __name__ == "__main__":
    admin_token = get_admin_access_token()
    user_map = create_users(admin_token)
    project_map = create_projects(admin_token)
    assign_project_roles(admin_token, user_map, project_map)

    print("Setup completed successfully!")
    print("\nTest Users:")
    for user_data in USERS_TO_CREATE:
        print(f"   {user_data['email']} (password: {user_data['password']})")

    print("\nTest Projects:")
    for project_name in project_map:
        print(f"   {project_name}")
