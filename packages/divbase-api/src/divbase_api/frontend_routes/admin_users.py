"""
Frontend routes for Admin to manage DivBase users.

All routes here should rely on get_current_admin_user_from_cookie dependency to ensure user is admin.
"""

import logging

from fastapi import APIRouter, Depends, Form, Request, status
from fastapi.responses import HTMLResponse, RedirectResponse
from pydantic import SecretStr
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.users import create_user, get_all_users, get_user_by_id
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
    current_admin: UserDB = Depends(get_current_admin_user_from_cookie),
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
            "current_user": UserResponse.model_validate(current_admin),
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
        context={"request": request, "current_user": UserResponse.model_validate(current_admin)},
    )


@fr_admin_users_router.post("/create", response_class=HTMLResponse)
async def post_create_user_endpoint(
    name: str = Form(...),
    email: str = Form(...),
    password: str = Form(...),
    is_admin: bool = Form(False),
    db: AsyncSession = Depends(get_db),
    current_admin: UserDB = Depends(get_current_admin_user_from_cookie),
):
    """Post request to create a new regular or admin user."""
    user_data = UserCreate(name=name, email=email, password=SecretStr(password))

    new_user = await create_user(db=db, user_data=user_data, is_admin=is_admin)
    logger.info(f"Admin user: {current_admin.email} created a new user: {new_user.email}")
    return RedirectResponse(url="/admin/users", status_code=status.HTTP_303_SEE_OTHER)


@fr_admin_users_router.get("/edit/{user_id}", response_class=HTMLResponse, status_code=status.HTTP_200_OK)
async def get_edit_user_endpoint(
    request: Request,
    user_id: int,
    current_admin: UserDB = Depends(get_current_admin_user_from_cookie),
    db: AsyncSession = Depends(get_db),
):
    """Page for editing a user."""
    user = await get_user_by_id(db=db, id=user_id)
    if not user:
        logger.warning(f"Admin user: {current_admin.email} tried to EDIT a non-existing user: {user_id}")
        return RedirectResponse(url="/admin/users", status_code=status.HTTP_303_SEE_OTHER)

    return templates.TemplateResponse(
        request=request,
        name="admin_pages/edit_user.html",
        context={
            "request": request,
            "current_user": UserResponse.model_validate(current_admin),
            "user_to_edit": UserResponse.model_validate(user),
        },
    )


# TODO post request of above.


@fr_admin_users_router.post("/{user_id}/inactivate", response_class=HTMLResponse)
async def post_inactivate_user_endpoint(
    user_id: int,
    db: AsyncSession = Depends(get_db),
    current_admin: UserDB = Depends(get_current_admin_user_from_cookie),
):
    """Post request to inactivate a user."""
    user = await get_user_by_id(db=db, id=user_id)
    if not user:
        logger.warning(f"Admin user: {current_admin.email} tried to DEACTIVATE a non-existing user: {user_id}")
        return RedirectResponse(url="/admin/users", status_code=status.HTTP_303_SEE_OTHER)

    user.is_active = False
    await db.commit()
    logger.info(f"Admin user: {current_admin.email} DEACTIVATED user: {user.email}")
    return RedirectResponse(url="/admin/users", status_code=status.HTTP_303_SEE_OTHER)


@fr_admin_users_router.post("/{user_id}/activate", response_class=HTMLResponse)
async def post_activate_user_endpoint(
    user_id: int,
    db: AsyncSession = Depends(get_db),
    current_admin: UserDB = Depends(get_current_admin_user_from_cookie),
):
    """Post request to activate a user."""
    user = await get_user_by_id(db=db, id=user_id)
    if not user:
        logger.warning(f"Admin user: {current_admin.email} tried to ACTIVATE a non-existing user: {user_id}")
        return RedirectResponse(url="/admin/users", status_code=status.HTTP_303_SEE_OTHER)

    user.is_active = True
    await db.commit()
    logger.info(f"Admin user: {current_admin.email} ACTIVATED user: {user.email}")
    return RedirectResponse(url="/admin/users", status_code=status.HTTP_303_SEE_OTHER)
