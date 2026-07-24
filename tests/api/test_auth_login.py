"""Auth API tests for login via password."""

from tests.api.helpers import TEST_PASSWORD, create_test_user

LOGIN_ENDPOINT = "/api/v1/auth/login"
GENERIC_LOGIN_FAILED_ERROR_MSG = "Invalid email or password or user account does not exist."


def _send_login_request(divbase_client, email: str, password: str = TEST_PASSWORD):
    return divbase_client.post(
        LOGIN_ENDPOINT,
        data={"username": email, "password": password},
    )


async def test_login_with_wrong_email_or_password(divbase_client, db_session):
    user = await create_test_user(db_session)
    response = await _send_login_request(divbase_client, email=user.email, password="wrong-password")
    assert response.status_code == 401

    response = await _send_login_request(divbase_client, email="some-email@example.com")
    assert response.status_code == 401


async def test_login_as_invalid_user_is_rejected(divbase_client, db_session):
    """3 diff ways to be classified as an invalid user: email not verified, deactivated, soft deleted."""
    deactivated_user = await create_test_user(db_session, is_active=False)
    soft_deleted_user = await create_test_user(db_session, is_soft_deleted=True)
    email_not_verified_user = await create_test_user(db_session, email_verified=False)

    response = await _send_login_request(divbase_client, email=deactivated_user.email)
    assert response.status_code == 401
    assert response.json()["detail"] == GENERIC_LOGIN_FAILED_ERROR_MSG

    response = await _send_login_request(divbase_client, email=soft_deleted_user.email)
    assert response.status_code == 401
    assert response.json()["detail"] == GENERIC_LOGIN_FAILED_ERROR_MSG

    # email_not_verified is a special case:
    # if user gave correct password, we tell them they need to verfiy their email
    # if wrong password, we give the generic error to avoid email enumeration
    response = await _send_login_request(divbase_client, email=email_not_verified_user.email)
    assert response.status_code == 401
    assert "not verified" in response.json()["detail"].lower()

    response = await _send_login_request(divbase_client, email=email_not_verified_user.email, password="wrongpassword")
    assert response.status_code == 401
    assert response.json()["detail"] == GENERIC_LOGIN_FAILED_ERROR_MSG
