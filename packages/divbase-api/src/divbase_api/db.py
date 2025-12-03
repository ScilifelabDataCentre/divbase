"""
Handles connection between FastAPI and the postgresql db.
"""

import logging
from pathlib import Path
from typing import AsyncGenerator

from alembic import command
from alembic.config import Config
from alembic.util.exc import CommandError
from sqlalchemy import text
from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker, create_async_engine

from divbase_api.api_config import settings
from divbase_api.crud.users import create_user, get_all_users
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


def check_db_migrations_up_to_date() -> None:
    """
    Checks if the DB is at HEAD (has applied all migrations) and raise an error if not.
    """
    alembic_config_path = Path(__file__).parent / "alembic.ini"
    try:
        alembic_cfg = Config(file_=alembic_config_path)
        command.check(alembic_cfg)

    except CommandError:
        logger.error(
            """
            It seems like database migrations are pending. You need to run them before continuing.

            - If you are doing local development:
            Restarting the stack will automatically apply the migrations.
            (Down and up the containers, compose watch does not cover this)
            migrations are run in a docker init container

            - If you are in production/deployed environments:
            WARNING: Running migrations in a deployed environment is a dangerous operation. 
            1. Scale down the API and worker replicas to 0. 
            2. Run the migration job (ensure you take a database backup first).
            3. Scale the services back up after the migration is complete.
            Exiting...
            """
        )
        raise


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
            email_verified=True,
        )
    logger.info(f"First admin user created with email: {admin_user.email}")
