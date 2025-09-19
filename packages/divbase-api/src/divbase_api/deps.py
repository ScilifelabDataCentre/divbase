"""
Authentication for FastAPI routes.

TODO - perhaps token types should be a str enum
"""

from typing import Annotated

from fastapi import Depends, HTTPException, status
from fastapi.security import OAuth2PasswordBearer
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.users import get_user_by_id
from divbase_api.db import get_db
from divbase_api.models.users import UserDB
from divbase_api.security import verify_token

oauth2_scheme = OAuth2PasswordBearer(tokenUrl="/api/v1/auth/login")


async def get_current_user(token: Annotated[str, Depends(oauth2_scheme)], db: AsyncSession = Depends(get_db)) -> UserDB:
    """
    Get current user from JWT Access token

    We should not be specific about why/if credentials are invalid.
    """
    credentials_exception = HTTPException(
        status_code=status.HTTP_401_UNAUTHORIZED,
        detail="Could not validate credentials",
        headers={"WWW-Authenticate": "Bearer"},
    )

    user_id = verify_token(token, desired_token_type="access")
    if user_id is None:
        raise credentials_exception

    user = await get_user_by_id(db=db, id=user_id)
    if user is None or not user.is_active:
        raise credentials_exception

    return user


async def get_current_admin_user(current_user: Annotated[UserDB, Depends(get_current_user)]) -> UserDB:
    """
    Verify current user has admin access.

    Note: This function works by using a sub-dependency to get the current user,
    so all the checks in the sub-dependancy function are ran when you use this dependency.
    (see here: https://fastapi.tiangolo.com/yo/advanced/security/oauth2-scopes/#dependency-tree-and-scopes)
    """
    if not current_user.is_admin:
        raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Admin access required")
    return current_user
