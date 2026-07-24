"""Auth API tests for refresh tokens"""

from datetime import datetime, timezone

import structlog
from pydantic import SecretStr

from divbase_api.crud.auth import update_user_password
from divbase_api.crud.users import get_user_by_id_or_raise
from divbase_api.schemas.users import UserPasswordUpdate
from divbase_lib.api_schemas.auth import RefreshTokenResponse
from tests.api.helpers import create_test_user, login


async def test_refresh_with_valid_token_returns_new_access_token(divbase_client, db_session):
    user = await create_test_user(db_session)
    login_data = await login(divbase_client, user.email)

    response = await divbase_client.post(
        "/api/v1/auth/refresh",
        json={"refresh_token": login_data.refresh_token},
    )
    assert response.status_code == 200
    refresh_data = RefreshTokenResponse.model_validate(response.json())
    assert refresh_data.access_token != login_data.access_token


async def test_refresh_with_garbage_token_is_rejected(divbase_client):
    response = await divbase_client.post(
        "/api/v1/auth/refresh",
        json={"refresh_token": "not-a-real-token"},
    )
    assert response.status_code == 401


async def test_refresh_with_access_token_instead_of_refresh_token_is_rejected(divbase_client, db_session):
    """The refresh endpoint must reject a token of the wrong type, even if otherwise validly signed."""
    user = await create_test_user(db_session)
    login_data = await login(divbase_client, user.email)

    response = await divbase_client.post(
        "/api/v1/auth/refresh",
        json={"refresh_token": login_data.access_token},
    )
    assert response.status_code == 401


async def test_refresh_token_is_revoked_after_logout(divbase_client, db_session):
    """Can't reuse a refresh token as should be revoked on logout."""
    user = await create_test_user(db_session)
    login_data = await login(divbase_client, user.email)

    logout_response = await divbase_client.post(
        "/api/v1/auth/logout",
        json={"refresh_token": login_data.refresh_token},
    )
    assert logout_response.status_code == 204

    with structlog.testing.capture_logs() as cap_logs:
        refresh_response = await divbase_client.post(
            "/api/v1/auth/refresh",
            json={"refresh_token": login_data.refresh_token},
        )
    assert refresh_response.status_code == 401
    assert any("Attempt to use revoked refresh token" in entry.get("event", "") for entry in cap_logs)


async def test_refresh_token_issued_before_password_change_is_rejected(divbase_client, db_session):
    """If a user changes their password, any refresh tokens issued before that change should be rejected."""
    user = await create_test_user(db_session)
    login_data = await login(divbase_client, user.email)

    new_password_data = UserPasswordUpdate(password=SecretStr("newpassword"), confirm_password=SecretStr("newpassword"))
    await update_user_password(db=db_session, user_id=user.id, password_data=new_password_data)

    response = await divbase_client.post(
        "/api/v1/auth/refresh",
        json={"refresh_token": login_data.refresh_token},
    )
    assert response.status_code == 401


async def test_refresh_token_belonging_to_an_invalid_user_is_rejected(divbase_client, db_session):
    """3 diff ways to be classified as an invalid user: email not verified, deactivated, soft deleted."""
    # In this case we need to login first to get a refresh token, then make the user invalid, otherwise login will fail.
    deactivated_user = await create_test_user(db_session)
    soft_deleted_user = await create_test_user(db_session)
    email_not_verified_user = await create_test_user(db_session)

    deactivated_login_data = await login(divbase_client, email=deactivated_user.email)
    soft_deleted_login_data = await login(divbase_client, email=soft_deleted_user.email)
    email_not_verified_login_data = await login(divbase_client, email=email_not_verified_user.email)

    deactivated_user_row = await get_user_by_id_or_raise(db=db_session, id=deactivated_user.id)
    deactivated_user_row.is_active = False

    soft_deleted_user_row = await get_user_by_id_or_raise(db=db_session, id=soft_deleted_user.id)
    soft_deleted_user_row.is_deleted = True
    soft_deleted_user_row.date_deleted = datetime.now(timezone.utc)

    email_not_verified_user_row = await get_user_by_id_or_raise(db=db_session, id=email_not_verified_user.id)
    email_not_verified_user_row.email_verified = False

    for login_data in [deactivated_login_data, soft_deleted_login_data, email_not_verified_login_data]:
        response = await divbase_client.post(
            "/api/v1/auth/refresh",
            json={"refresh_token": login_data.refresh_token},
        )
        assert response.status_code == 401
