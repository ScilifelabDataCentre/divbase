"""
Fixtures for API tests.

The fixtures build the FastAPI app with a test db instance, and provide an async test client for sending requests to the app.
"""

from pathlib import Path

import pytest
import pytest_asyncio
from alembic import command
from alembic.config import Config
from celery.backends.database import DatabaseBackend
from httpx import ASGITransport, AsyncClient
from pytest_mock_resources import create_postgres_fixture
from sqlalchemy.ext.asyncio import AsyncSession, create_async_engine

from divbase_api.db import get_db
from divbase_api.divbase_api import create_app
from divbase_api.worker.tasks import app as celery_app

ALEMBIC_CONFIG_PATH = (
    Path(__file__).parent.parent.parent / "packages" / "divbase-api" / "src" / "divbase_api" / "alembic.ini"
)


def _migrate_db_to_head(conn) -> None:
    """Run alembic migrations + add the Celery managed tables on the test postgres db."""
    alembic_cfg = Config(file_=str(ALEMBIC_CONFIG_PATH))
    alembic_cfg.set_main_option("sqlalchemy.url", conn.engine.url.render_as_string(hide_password=False))
    command.upgrade(alembic_cfg, "head")

    # Celery manages these tables, not alembic, so they are not created by the migrations above.)
    DatabaseBackend(app=celery_app, url=conn.engine.url.render_as_string(hide_password=False))


test_pg_engine = create_postgres_fixture(_migrate_db_to_head, scope="session")


@pytest_asyncio.fixture()
async def db_session(test_pg_engine):
    """DB session for API tests, automatically rolls back after each test."""
    async_url = test_pg_engine.url.render_as_string(hide_password=False).replace(
        "postgresql+psycopg2", "postgresql+asyncpg"
    )
    async_engine = create_async_engine(async_url)
    async with async_engine.connect() as conn:
        await conn.begin()
        session = AsyncSession(bind=conn, join_transaction_mode="create_savepoint", expire_on_commit=False)
        try:
            yield session
        finally:
            await session.close()
            await conn.rollback()
    await async_engine.dispose()


@pytest.fixture
def app(db_session):
    """DivBase FastAPI app without lifespan startup checks, using the test db_session."""
    fastapi_app = create_app(lifespan=None)
    fastapi_app.dependency_overrides[get_db] = lambda: db_session
    return fastapi_app


@pytest_asyncio.fixture
async def divbase_client(app):
    """Async HTTP client that tests can use to send requests to divbase-api"""
    transport = ASGITransport(app=app)
    async with AsyncClient(transport=transport, base_url="http://localhost") as divbase_client:
        yield divbase_client
