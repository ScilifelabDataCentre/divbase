"""
Handles connection between FastAPI and the postgresql db.
"""

import logging
from typing import AsyncGenerator

from sqlalchemy import text
from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker, create_async_engine

from divbase_api.api_config import settings
from divbase_api.crud.users import create_user, get_all_users
from divbase_api.models.base import Base
from divbase_api.schemas.users import UserCreate

logger = logging.getLogger(__name__)


# NOTE: The creation of the AsyncSessionLocal should only occur once, when the module is loaded.
# AKA: Do not refactor to have these in functions.
engine = create_async_engine(
    url=settings.database.url.get_secret_value(),
    echo=settings.database.echo_db_output,
    pool_pre_ping=True,
)
AsyncSessionLocal = async_sessionmaker(
    engine,
    class_=AsyncSession,
    expire_on_commit=False,
)


async def get_db() -> AsyncGenerator[AsyncSession, None]:
    """Dependency to get database session. To be used in FastAPI endpoints."""
    async with AsyncSessionLocal() as session:
        try:
            yield session
        finally:
            await session.close()


async def health_check_db() -> bool:
    """Check if the database connection is healthy."""
    try:
        async with engine.begin() as conn:
            await conn.execute(text("SELECT 1"))
        return True
    except Exception:
        return False


async def create_all_tables() -> None:
    """Create all tables in the db."""
    async with engine.begin() as conn:
        await conn.run_sync(Base.metadata.create_all)


async def drop_all_tables() -> None:
    """Drop all tables in the db."""
    async with engine.begin() as conn:
        await conn.run_sync(Base.metadata.drop_all)


async def create_first_admin_user() -> None:
    """
    Create the first admin user if:
    1. no admin users already exists in the db.
    2. FIRST_ADMIN_EMAIL and FIRST_ADMIN_PASSWORD env vars are set.
    """
    admin_email = settings.api.first_admin_email
    admin_password = settings.api.first_admin_password

    if admin_email == "NOT_SET" or admin_password.get_secret_value() == "NOT_SET":
        logger.info(
            "No first admin env vars set (FIRST_ADMIN_EMAIL, FIRST_ADMIN_PASSWORD), skipping first admin creation"
        )
        return

    async with AsyncSessionLocal() as db:
        admin_users = await get_all_users(db=db, admins_only=True)

        if len(admin_users) > 0:
            logger.info("At least 1 admin user already exists, skipping first admin user auto creation")
            return

        user_info = UserCreate(
            name="First Admin",
            email=admin_email,
            password=admin_password,
        )

        admin_user = await create_user(
            db=db,
            user_data=user_info,
            is_admin=True,
        )
    logger.info(f"First admin user created with email: {admin_user.email}")
