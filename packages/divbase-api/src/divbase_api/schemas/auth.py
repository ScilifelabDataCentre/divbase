"""
Schemas for login + access and refresh tokens
"""

from pydantic import BaseModel, Field


class CLILoginResponse(BaseModel):
    """Response model for API (aka divbase-cli) login endpoint."""

    access_token: str = Field(..., description="Bearer access token for authentication")
    access_token_expires_at: int = Field(..., description="Unix timestamp when the access token expires")
    refresh_token: str = Field(..., description="Bearer refresh token for obtaining new access tokens")
    refresh_token_expires_at: int = Field(..., description="Unix timestamp when the refresh token expires")
    email: str = Field(..., description="Email of the authenticated user")


class RefreshTokenRequest(BaseModel):
    """Request model for refresh token endpoint."""

    refresh_token: str = Field(..., description="Bearer refresh token for obtaining a new access token")


class RefreshTokenResponse(BaseModel):
    """Response model for refresh token endpoint."""

    access_token: str = Field(..., description="Bearer access token for authentication")
    expires_at: int = Field(..., description="Unix timestamp when the access token expires")


class LogoutRequest(BaseModel):
    """Request model for logout endpoint."""

    refresh_token: str = Field(..., description="Bearer refresh token to be revoked on logout")
