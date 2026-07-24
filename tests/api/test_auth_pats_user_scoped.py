"""Tests using Personal Access Tokens (PATs) for user scoped operations to DivBase API."""

from datetime import datetime, timedelta, timezone

from divbase_api.crud.personal_access_tokens import get_users_personal_access_tokens, soft_delete_personal_access_token
from divbase_api.schemas.users import UserResponse
from tests.api.helpers import bearer_header, create_pat, create_test_user

WHO_AMI_ENDPOINT = "/api/v1/auth/whoami"


async def test_whoami_with_valid_pat_returns_current_user(divbase_client, db_session):
    user = await create_test_user(db_session)
    raw_pat = await create_pat(db_session, user_id=user.id)

    response = await divbase_client.get(WHO_AMI_ENDPOINT, headers=bearer_header(raw_pat))

    assert response.status_code == 200
    whoami_data = UserResponse.model_validate(response.json())
    assert whoami_data.email == user.email


async def test_expired_or_deleted_pat_is_rejected(divbase_client, db_session):
    user = await create_test_user(db_session)

    deleted_pat = await create_pat(db_session, user_id=user.id)
    pats = await get_users_personal_access_tokens(db_session, user_id=user.id)
    await soft_delete_personal_access_token(db_session, pat_id=pats[0].id, user_id=user.id)
    deleted_response = await divbase_client.get(WHO_AMI_ENDPOINT, headers=bearer_header(deleted_pat))
    assert deleted_response.status_code == 401

    expired_pat = await create_pat(
        db_session, user_id=user.id, expires_at=datetime.now(timezone.utc) - timedelta(days=1)
    )
    expired_response = await divbase_client.get(WHO_AMI_ENDPOINT, headers=bearer_header(expired_pat))
    assert expired_response.status_code == 401


async def test_pat_belonging_to_an_invalid_user_is_rejected(divbase_client, db_session):
    """3 diff ways to be classified as an invalid user: email not verified, deactivated, soft deleted."""
    deactivated_user = await create_test_user(db_session, is_active=False)
    soft_deleted_user = await create_test_user(db_session, is_soft_deleted=True)
    email_not_verified_user = await create_test_user(db_session, email_verified=False)

    deactivated_pat = await create_pat(db_session, user_id=deactivated_user.id)
    soft_deleted_pat = await create_pat(db_session, user_id=soft_deleted_user.id)
    email_not_verified_pat = await create_pat(db_session, user_id=email_not_verified_user.id)

    for raw_pat in [deactivated_pat, soft_deleted_pat, email_not_verified_pat]:
        response = await divbase_client.get(WHO_AMI_ENDPOINT, headers=bearer_header(raw_pat))
        assert response.status_code == 401


async def test_task_history_endpoint_denied_for_pat_without_task_history_scope(divbase_client, db_session):
    user = await create_test_user(db_session)
    no_history_pat = await create_pat(db_session, user_id=user.id, permissions={"task_history": False})
    yes_history_pat = await create_pat(db_session, user_id=user.id, permissions={"task_history": True})

    response = await divbase_client.get("/api/v1/task-history/tasks/user", headers=bearer_header(no_history_pat))
    assert response.status_code == 403

    response = await divbase_client.get("/api/v1/task-history/tasks/user", headers=bearer_header(yes_history_pat))
    assert response.status_code == 200
