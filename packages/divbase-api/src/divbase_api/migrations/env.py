import os

from alembic import context
from sqlalchemy import engine_from_config, pool

from divbase_api.models.base import Base
from divbase_api.worker.tasks import CELERY_GROUPMETA_TABLE_NAME, CELERY_TASKMETA_TABLE_NAME

# this is the Alembic Config object, which provides
# access to the values within the .ini file in use.
config = context.config

# Interpret the config file for Python logging.
# This line sets up loggers basically.
# from logging.config import fileConfig
# if config.config_file_name is not None:
#     fileConfig(config.config_file_name, disable_existing_loggers=False)
# Commented out to use existing fastapi logging config

# add your model's MetaData object here
# Import your models' metadata
target_metadata = Base.metadata

# Set sqlalchemy.url
# Using the sync URL over async for simplicity (and no performance impact for something like this).
database_url = os.getenv("SYNC_DATABASE_URL")
if database_url:
    config.set_main_option("sqlalchemy.url", database_url)


def include_object(object, name, type_, reflected, compare_to) -> bool:
    """Exclude celery managed tables from auto generated migrations."""
    celery_managed_tables = {CELERY_TASKMETA_TABLE_NAME, CELERY_GROUPMETA_TABLE_NAME}
    return not (type_ == "table" and name in celery_managed_tables)


def run_migrations_offline() -> None:
    """Run migrations in 'offline' mode.

    This configures the context with just a URL
    and not an Engine, though an Engine is acceptable
    here as well.  By skipping the Engine creation
    we don't even need a DBAPI to be available.

    Calls to context.execute() here emit the given string to the
    script output.

    """
    url = config.get_main_option("sqlalchemy.url")
    context.configure(
        url=url,
        target_metadata=target_metadata,
        literal_binds=True,
        dialect_opts={"paramstyle": "named"},
        include_object=include_object,
    )

    with context.begin_transaction():
        context.run_migrations()


def run_migrations_online() -> None:
    """Run migrations in 'online' mode.

    In this scenario we need to create an Engine
    and associate a connection with the context.
    """
    connectable = context.config.attributes.get("connection", None)

    if connectable is None:
        connectable = engine_from_config(
            config.get_section(config.config_ini_section, {}),
            prefix="sqlalchemy.",
            poolclass=pool.NullPool,
        )

    with connectable.connect() as connection:
        context.configure(connection=connection, target_metadata=target_metadata, include_object=include_object)
        with context.begin_transaction():
            context.run_migrations()


if context.is_offline_mode():
    run_migrations_offline()
else:
    run_migrations_online()
