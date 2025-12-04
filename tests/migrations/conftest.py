"""
Pytest alembic used for automated migration testing.
https://github.com/schireson/pytest-alembic
"""

from pathlib import Path

import pytest
from pytest_alembic.config import Config
from pytest_mock_resources import create_postgres_fixture

# Set up pytest-alembic fixture for migrations
pg = create_postgres_fixture()


@pytest.fixture
def alembic_engine(pg):
    return pg


@pytest.fixture
def alembic_config():
    """Provide the Alembic configuration for pytest-alembic."""
    api_src_dir = Path(__file__).parent.parent.parent / "packages" / "divbase-api" / "src" / "divbase_api"
    ini_file_path = api_src_dir / "alembic.ini"
    migrations_dir_path = api_src_dir / "migrations"
    return Config(
        config_options={
            "file": str(ini_file_path),
            "script_location": str(migrations_dir_path),
        },
    )
