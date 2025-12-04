"""
Unit tests for the security module.
"""

from datetime import timedelta
from unittest.mock import patch

import pytest
from pydantic import SecretStr

from divbase_api.security import (
    TokenType,
    create_token,
    get_password_hash,
    token_expires_delta,
    verify_password,
    verify_token,
)

USER_ID = 1


@pytest.fixture(scope="module")
def valid_tokens():
    """
    Provide a dict of JWTs for testing purposes.

    "sub" of JWT is user_id
    """
    valid_tokens = {}
    valid_tokens[TokenType.ACCESS] = create_token(subject=USER_ID, token_type=TokenType.ACCESS).token
    valid_tokens[TokenType.REFRESH] = create_token(subject=USER_ID, token_type=TokenType.REFRESH).token
    valid_tokens[TokenType.EMAIL_VERIFICATION] = create_token(
        subject=USER_ID, token_type=TokenType.EMAIL_VERIFICATION
    ).token
    valid_tokens[TokenType.PASSWORD_RESET] = create_token(subject=USER_ID, token_type=TokenType.PASSWORD_RESET).token
    return valid_tokens


@pytest.fixture(scope="module")
def expired_tokens():
    expired_tokens = {}
    with patch.dict(
        token_expires_delta,
        {
            TokenType.ACCESS: timedelta(seconds=-10),
            TokenType.REFRESH: timedelta(seconds=-10),
            TokenType.EMAIL_VERIFICATION: timedelta(seconds=-10),
            TokenType.PASSWORD_RESET: timedelta(seconds=-10),
        },
    ):
        expired_tokens[TokenType.ACCESS] = create_token(subject=USER_ID, token_type=TokenType.ACCESS).token
        expired_tokens[TokenType.REFRESH] = create_token(subject=USER_ID, token_type=TokenType.REFRESH).token
        expired_tokens[TokenType.EMAIL_VERIFICATION] = create_token(
            subject=USER_ID, token_type=TokenType.EMAIL_VERIFICATION
        ).token
        expired_tokens[TokenType.PASSWORD_RESET] = create_token(
            subject=USER_ID, token_type=TokenType.PASSWORD_RESET
        ).token

    return expired_tokens


def test_validate_tokens(valid_tokens):
    """
    Test valid tokens can be validated and only validated by their correct type (access, refresh, etc).
    """
    for token_type, token in valid_tokens.items():
        token_data = verify_token(token, token_type)
        assert token_data is not None
        assert token_data.user_id == USER_ID


def test_expired_tokens_fail(expired_tokens):
    """
    Test that expired tokens are correctly identified as expired.
    """
    for token_type, token in expired_tokens.items():
        assert verify_token(token, token_type) is None


def test_invalid_tokens_fail():
    """Test that invalid tokens are correctly identified as invalid."""
    invalid_tokens = [
        "this.is.not.a.valid.token",
        "another.invalid.token",
        "",
        "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJzdWIiOjEsImlhdCI6MTY5ODUwMDAwMCwiZXhwIjoxNjk4NTAw",
    ]
    for token in invalid_tokens:
        for token_type in TokenType:
            assert verify_token(token=token, desired_token_type=token_type) is None


def test_incorrect_token_type_fails():
    """Test that tokens created for one type fail when attempting to validate as a different type."""
    for token_type in TokenType:
        token_data = create_token(subject=1, token_type=token_type)
        for test_type in TokenType:
            if test_type == token_type:
                token_data_verified = verify_token(token=token_data.token, desired_token_type=test_type)
                assert token_data_verified is not None
                assert token_data_verified.user_id == USER_ID
            else:
                assert verify_token(token=token_data.token, desired_token_type=test_type) is None


def test_password_hashing():
    """Test password hashing and verification."""
    password = "badpassword"
    hashed_password = get_password_hash(SecretStr(password))

    assert verify_password(plain_password=password, hashed_password=hashed_password)

    assert verify_password(plain_password="wrong_password", hashed_password=hashed_password) is False

    hashed_password_2 = get_password_hash(SecretStr(password))
    assert hashed_password != hashed_password_2
