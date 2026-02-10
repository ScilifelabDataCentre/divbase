"""
Pydantic schemas for user-related operations.

TODOs:
- Add proper password strength validation here.
"""

from datetime import datetime

from pydantic import (
    BaseModel,
    ConfigDict,
    EmailStr,
    Field,
    SecretStr,
    field_validator,
    model_validator,
)
from typing_extensions import Self


class UserBase(BaseModel):
    name: str = Field(..., min_length=3, max_length=100)
    email: EmailStr = Field(..., max_length=50)
    organisation: str = Field(..., max_length=200)
    organisation_role: str = Field(..., max_length=100)

    model_config = ConfigDict(from_attributes=True, str_strip_whitespace=True)

    @field_validator("email")
    @classmethod
    def lowercase_email(cls, v):
        return v.lower() if v else v


class UserCreate(UserBase):
    password: SecretStr = Field(..., min_length=8, max_length=100)


class UserUpdate(BaseModel):
    """Schema for a user to update their own profile."""

    name: str = Field(..., min_length=3, max_length=100)
    organisation: str = Field(..., max_length=200)
    organisation_role: str = Field(..., max_length=100)


class UserPasswordUpdate(BaseModel):
    """Schema for a user to update their own password."""

    password: SecretStr
    confirm_password: SecretStr

    @field_validator("password")
    @classmethod
    def password_strength(cls, v):
        if len(v.get_secret_value()) < 8:
            raise ValueError("Password must be at least 8 characters long")
        return v

    @model_validator(mode="after")
    def check_passwords_match(self) -> Self:
        if self.password.get_secret_value() != self.confirm_password.get_secret_value():
            raise ValueError("Passwords do not match")
        return self


class UserResponse(UserBase):
    """Response schema, aka returned by API endpoints."""

    id: int
    is_admin: bool
    is_active: bool
    created_at: datetime
    updated_at: datetime
    is_deleted: bool
