"""
The API server for DivBase.
"""

import logging
import sys
from contextlib import asynccontextmanager
from pathlib import Path
from typing import AsyncGenerator

import uvicorn
from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles

from divbase_api import __version__ as divbase_version
from divbase_api.admin_panel import register_admin_panel
from divbase_api.api_config import LOCAL_DEV_ENVIRONMENTS, api_settings
from divbase_api.db import (
    check_db_migrations_up_to_date,
    create_first_admin_user,
    engine,
    health_check_db,
)
from divbase_api.exception_handlers import register_exception_handlers
from divbase_api.frontend_routes.auth import fr_auth_router
from divbase_api.frontend_routes.core import fr_core_router
from divbase_api.frontend_routes.profile import fr_profile_router
from divbase_api.frontend_routes.projects import fr_projects_router
from divbase_api.middleware import register_middleware
from divbase_api.routes.admin import admin_router
from divbase_api.routes.auth import auth_router
from divbase_api.routes.core import core_router
from divbase_api.routes.project_versions import project_version_router
from divbase_api.routes.queries import query_router
from divbase_api.routes.s3 import s3_router
from divbase_api.routes.task_history import task_history_router
from divbase_api.routes.vcf_dimensions import vcf_dimensions_router

logging.basicConfig(level=api_settings.general.log_level, handlers=[logging.StreamHandler(sys.stderr)])

logger = logging.getLogger(__name__)


@asynccontextmanager
async def lifespan(app: FastAPI) -> AsyncGenerator[None, None]:
    """Lifespan context manager that handles startup and shutdown of API server"""
    # startup
    logger.info("Starting up DivBase API...")

    api_settings.validate_api_settings()
    logger.info("All API settings are correctly set.")

    if not await health_check_db():
        raise ConnectionError("Could not connect to the database or db unhealthy. Exiting...")
    logger.info("Database connection healthy.")

    check_db_migrations_up_to_date()
    logger.info("db is up to date with the latest migration files.")

    await create_first_admin_user()
    logger.info("DivBase API startup events complete.")

    yield

    # cleanup
    logger.info("Shutting down, closing any DB connections")
    await engine.dispose()


app = FastAPI(
    lifespan=lifespan,
    title="DivBase API",
    docs_url="/api/v1/docs",
    version=divbase_version,
    summary="DivBase API is used by divbase-cli to interact with the DivBase. We do not recommend users to interact directly with the API.",
    description="""
    Users are strongly encouraged to use divbase-cli rather than calling the API endpoints directly
    
    divbase-cli offers several advantages over direct API usage, including:
    - Handles authentication logic, including automatic token refresh
    - Simplifies commands and workflows (e.g. file uploads/downloads require working with presigned URLs).
    - Better error messages and user experience

    If there is something you cannot do with divbase-cli that you think should be possible, please let us know.

    Visit our docs site for more info on how to use divbase-cli: https://scilifelabdatacentre.github.io/divbase/
    """,
)

app.include_router(auth_router, prefix="/api/v1/auth", tags=["auth"])
app.include_router(core_router, prefix="/api/v1/core", tags=["core"])
app.include_router(project_version_router, prefix="/api/v1/project-versions", tags=["project-versioning"])
app.include_router(query_router, prefix="/api/v1/query", tags=["query"])
app.include_router(s3_router, prefix="/api/v1/s3", tags=["s3"])
app.include_router(task_history_router, prefix="/api/v1/task-history", tags=["task-history"])
app.include_router(vcf_dimensions_router, prefix="/api/v1/vcf-dimensions", tags=["vcf-dimensions"])
if api_settings.general.environment in LOCAL_DEV_ENVIRONMENTS:
    # not needed in deployed enviroments, so no need to expose it.
    app.include_router(admin_router, prefix="/api/v1/admin", tags=["admin"])

app.include_router(fr_auth_router, prefix="", include_in_schema=False)
app.include_router(fr_core_router, prefix="", include_in_schema=False)
app.include_router(fr_profile_router, prefix="/profile", include_in_schema=False)
app.include_router(fr_projects_router, prefix="/projects", include_in_schema=False)


register_exception_handlers(app)
register_admin_panel(app=app, engine=engine)
register_middleware(app=app)

static_dir_path = Path(__file__).parent / "static"
app.mount("/static", StaticFiles(directory=static_dir_path), name="static")


def main():
    uvicorn.run("divbase_api.divbase_api:app", host="127.0.0.1", port=8000, reload=True)


if __name__ == "__main__":
    main()
