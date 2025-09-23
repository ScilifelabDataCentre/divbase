"""
Pydantic schemas for user-related operations.

TODOs:
- Add proper password strength validation here.
"""

from datetime import datetime

from pydantic import BaseModel, ConfigDict, EmailStr, Field, SecretStr


class UserBase(BaseModel):
    name: str = Field(..., min_length=3, max_length=100)
    email: EmailStr

    model_config = ConfigDict(from_attributes=True, str_strip_whitespace=True)


class UserCreate(UserBase):
    password: SecretStr = Field(..., min_length=8, max_length=100)


class UserUpdate(BaseModel):
    name: str | None = Field(None, min_length=3, max_length=100)
    email: EmailStr | None = None


class UserPasswordUpdate(BaseModel):
    password: SecretStr | None = Field(None, min_length=8, max_length=100)
    confirm_password: SecretStr | None = Field(None, min_length=8, max_length=100)


class UserResponse(UserBase):
    """Response schema, aka returned by API endpoints."""

    id: int
    is_admin: bool
    is_active: bool
    created_at: datetime
    updated_at: datetime
