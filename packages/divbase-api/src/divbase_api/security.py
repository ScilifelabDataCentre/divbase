"""
Security utilities for handling passwords and (TODO) JWT tokens.

bcrypt recommended by FastAPI docs: https://fastapi.tiangolo.com/yo/tutorial/security/oauth2-jwt/
Follows setup from official full stack template: https://github.com/fastapi/full-stack-fastapi-template/blob/master/backend/app/core/security.py
"""

from passlib.context import CryptContext
from pydantic import SecretStr

pwd_context = CryptContext(schemes=["bcrypt"], deprecated="auto")


def verify_password(plain_password: str, hashed_password: str) -> bool:
    return pwd_context.verify(plain_password, hashed_password)


def get_password_hash(password: SecretStr) -> str:
    return pwd_context.hash(password.get_secret_value())
