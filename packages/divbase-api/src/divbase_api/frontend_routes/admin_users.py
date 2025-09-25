"""
Frontend routes for Admin to manage DivBase users.

All routes here should rely on get_current_admin_user_from_cookie dependency to ensure user is admin.
"""

import logging

from fastapi import APIRouter, Depends, Form, Query, Request, status
from fastapi.responses import HTMLResponse, RedirectResponse
from pydantic import SecretStr
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.users import create_user, get_all_users
from divbase_api.db import get_db
from divbase_api.deps import get_current_admin_user_from_cookie
from divbase_api.frontend_routes.core import templates
from divbase_api.models.users import UserDB
from divbase_api.schemas.users import UserCreate, UserResponse

logger = logging.getLogger(__name__)

fr_admin_users_router = APIRouter()


@fr_admin_users_router.get("/", response_class=HTMLResponse)
async def admin_users_manage_endpoint(
    request: Request,
    current_admin_user: UserDB = Depends(get_current_admin_user_from_cookie),
    db: AsyncSession = Depends(get_db),
):
    """
    Render the admin's users management page.

    TODO - pagination
    """
    all_users = await get_all_users(db=db, limit=100)

    return templates.TemplateResponse(
        request=request,
        name="admin_pages/users.html",
        context={
            "request": request,
            "user": UserResponse.model_validate(current_admin_user),
            "users": [UserResponse.model_validate(user) for user in all_users],
        },
    )


@fr_admin_users_router.get("/create", response_class=HTMLResponse, status_code=status.HTTP_200_OK)
async def get_create_user_endpoint(
    request: Request,
    current_admin: UserDB = Depends(get_current_admin_user_from_cookie),
):
    """Page for creating a new regular or admin user."""
    return templates.TemplateResponse(
        request=request,
        name="admin_pages/create_user.html",
        context={"request": request, "user": UserResponse.model_validate(current_admin)},
    )


@fr_admin_users_router.post("/create", response_model=UserResponse, status_code=status.HTTP_303_SEE_OTHER)
async def post_create_user_endpoint(
    name: str = Form(...),
    email: str = Form(...),
    password: str = Form(...),
    is_admin: bool = Query(False, description="Set to true to create an admin user"),
    db: AsyncSession = Depends(get_db),
    current_admin: UserDB = Depends(get_current_admin_user_from_cookie),
):
    """Post request to create a new regular or admin user."""
    user_data = UserCreate(name=name, email=email, password=SecretStr(password))

    new_user = await create_user(db=db, user_data=user_data, is_admin=is_admin)
    logger.info(f"Admin user: {current_admin.email} created a new user: {new_user.email}")
    return RedirectResponse(url="/admin/users")
