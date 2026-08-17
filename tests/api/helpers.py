"""Shared helpers for DivBase API tests: e.g.: create users, projects, memberships, personal access tokens and logging in."""

from datetime import datetime, timedelta, timezone
from uuid import uuid4

from httpx import AsyncClient
from pydantic import SecretStr
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.personal_access_tokens import create_personal_access_token
from divbase_api.crud.projects import add_project_member, create_project
from divbase_api.crud.users import create_user, get_user_by_id_or_raise
from divbase_api.models.projects import ProjectRoles
from divbase_api.schemas.projects import ProjectCreate, ProjectResponse
from divbase_api.schemas.users import UserCreate, UserResponse
from divbase_lib.api_schemas.auth import CLILoginResponse

TEST_PASSWORD = "badpassword"


async def create_test_user(
    db_session: AsyncSession,
    email: str | None = None,
    password: str = TEST_PASSWORD,
    email_verified: bool = True,
    is_admin: bool = False,
    is_active: bool = True,
    is_soft_deleted: bool = False,
) -> UserResponse:
    """
    Creates a user for testing.

    Last password change is set to be one hour ago. This is because
    verify_user_from_refresh_token() rejects any refresh token issued before last password change (seconds precision).
    """
    if email is None:
        email = f"{uuid4()}@divbase.com"

    user_data = UserCreate(
        name="Test User",
        email=email,
        organisation="DivBase University",
        organisation_role="Developer",
        password=SecretStr(password),
        confirm_password=SecretStr(password),
    )
    user = await create_user(db=db_session, user_data=user_data, is_admin=is_admin, email_verified=email_verified)

    user_row = await get_user_by_id_or_raise(db=db_session, id=user.id)
    user_row.last_password_change = datetime.now(timezone.utc) - timedelta(hours=1)
    if not is_active:
        user_row.is_active = False
    if is_soft_deleted:
        user_row.is_deleted = True
        user_row.date_deleted = datetime.now(timezone.utc)

    await db_session.commit()
    return user


def bearer_header(token: str) -> dict[str, str]:
    """Create an authorization header for a JWT access or refresh token or a Personal Access Token (PAT)."""
    return {"Authorization": f"Bearer {token}"}


async def login(divbase_client: AsyncClient, email: str, password: str = TEST_PASSWORD) -> CLILoginResponse:
    """Log in to DivBase via the API endpoint, returns the CLILoginResponse"""
    response = await divbase_client.post(
        "/api/v1/auth/login",
        data={"username": email, "password": password},
    )
    assert response.status_code == 200
    return CLILoginResponse.model_validate(response.json())


async def create_pat(
    db_session: AsyncSession,
    user_id: int,
    name: str | None = None,
    permissions: dict | None = None,
    expires_at: datetime | None = None,
) -> str:
    """Create a PAT for a user, returning the token to use as a Bearer token in API requests."""
    if not name:
        name = f"test-pat-{uuid4()}"

    _, raw_pat = await create_personal_access_token(
        db=db_session,
        user_id=user_id,
        name=name,
        permissions=permissions or {},
        expires_at=expires_at,
    )
    return raw_pat.get_secret_value()


async def create_test_project(
    db_session,
    name: str = "test-project",
    storage_quota_bytes: int = 1_000_000,
) -> ProjectResponse:
    proj_data = ProjectCreate(
        name=name,
        description="A test project",
        bucket_name=f"{name}-bucket",
        storage_quota_bytes=storage_quota_bytes,
    )
    return await create_project(db=db_session, proj_data=proj_data)


async def add_test_project_member(db_session, project_id: int, user_email: str, role: ProjectRoles) -> None:
    await add_project_member(db=db_session, project_id=project_id, user_email=user_email, role=role)
