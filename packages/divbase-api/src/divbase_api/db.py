"""
Handles connection between FastAPI and the postgresql db.

# TODO:
- Add db health check
- Add ability to create first admin user if none exist.
- Read into async_sessionmaker params more e.g. pool_recycle.
"""

from typing import AsyncGenerator

from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker, create_async_engine

from divbase_api.config import settings
from divbase_api.models.base import Base

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


async def create_all_tables() -> None:
    """Create all tables in the db."""
    async with engine.begin() as conn:
        await conn.run_sync(Base.metadata.create_all)


async def drop_all_tables() -> None:
    """Drop all tables in the db."""
    async with engine.begin() as conn:
        await conn.run_sync(Base.metadata.drop_all)
