"""
Pydantic schemas for user-related operations.

TODOs:
- Add proper password strength validation here.
"""

from datetime import datetime

from pydantic import BaseModel, ConfigDict, EmailStr, Field, SecretStr, field_validator


class UserBase(BaseModel):
    name: str = Field(..., min_length=3, max_length=100)
    email: EmailStr = Field(..., max_length=50)

    model_config = ConfigDict(from_attributes=True, str_strip_whitespace=True)

    @field_validator("email")
    @classmethod
    def lowercase_email(cls, v):
        return v.lower() if v else v


class UserCreate(UserBase):
    password: SecretStr = Field(..., min_length=8, max_length=100)


class UserUpdate(UserBase):
    """Schema for a user to update their own profile."""

    pass


class AdminUserUpdate(UserBase):
    """Schema for admins to update a user."""

    is_admin: bool | None = None
    is_active: bool | None = None


class UserPasswordUpdate(BaseModel):
    """Schema for a user to update their own password."""

    current_password: SecretStr = Field(..., min_length=8, max_length=100)
    password: SecretStr | None = Field(None, min_length=8, max_length=100)
    confirm_password: SecretStr | None = Field(None, min_length=8, max_length=100)


class AdminUserPasswordUpdate(BaseModel):
    """Schema for an admin to update a user's password."""

    password: SecretStr | None = Field(None, min_length=8, max_length=100)
    confirm_password: SecretStr | None = Field(None, min_length=8, max_length=100)


class UserResponse(UserBase):
    """Response schema, aka returned by API endpoints."""

    id: int
    is_admin: bool
    is_active: bool
    created_at: datetime
    updated_at: datetime
