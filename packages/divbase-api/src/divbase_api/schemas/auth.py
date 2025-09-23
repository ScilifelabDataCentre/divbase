"""
Schemas for login + access and refresh tokens
"""

from pydantic import BaseModel, Field


class CLILoginResponse(BaseModel):
    """Response model for API (aka divbase-cli) login endpoint."""

    access_token: str = Field(..., description="Bearer access token for authentication")
    refresh_token: str = Field(..., description="Bearer refresh token for obtaining new access tokens")
    email: str = Field(..., description="Email of the authenticated user")


class WebLoginResponse(BaseModel):
    """Response model for web login endpoint."""

    pass  # TODO


class RefreshTokenRequest(BaseModel):
    """Request model for refresh token endpoint."""

    refresh_token: str = Field(..., description="Bearer refresh token for obtaining a new access token")


class RefreshTokenResponse(BaseModel):
    """Response model for refresh token endpoint."""

    access_token: str = Field(..., description="Bearer access token for authentication")
