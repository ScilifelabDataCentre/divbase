"""
The API server for DivBase.

TODO: user_name would later be determined by the authentication system.
"""

import logging
import sys
from contextlib import asynccontextmanager
from typing import AsyncGenerator

import uvicorn
from fastapi import FastAPI

from divbase_api.config import settings
from divbase_api.db import (
    create_all_tables,
    create_first_admin_user,
    engine,
    health_check_db,
)
from divbase_api.exception_handlers import register_exception_handlers
from divbase_api.frontend_routes.admin import fr_admin_router
from divbase_api.frontend_routes.auth import fr_auth_router
from divbase_api.frontend_routes.core import fr_core_router
from divbase_api.frontend_routes.profile import fr_profile_router
from divbase_api.get_task_history import get_task_history
from divbase_api.routes.admin import admin_router
from divbase_api.routes.auth import auth_router
from divbase_api.routes.projects import projects_router
from divbase_api.routes.users import users_router
from divbase_worker.tasks import (
    bcftools_pipe_task,
    sample_metadata_query_task,
    update_vcf_dimensions_task,
)

logging.basicConfig(level=settings.api.log_level, handlers=[logging.StreamHandler(sys.stdout)])

logger = logging.getLogger(__name__)


@asynccontextmanager
async def lifespan(app: FastAPI) -> AsyncGenerator[None, None]:
    """Lifespan context manager that handles startup and shutdown of API server"""
    # startup
    logger.info("Starting up DivBase API...")

    if not await health_check_db():
        raise ConnectionError("Could not connect to the database or db unhealthy. Exiting...")
    logger.info("Database connection healthy.")
    await create_all_tables()
    await create_first_admin_user()
    logger.info("DivBase API startup events complete.")

    yield

    # cleanup
    logger.info("Shutting down, closing any DB connections")
    await engine.dispose()


app = FastAPI(lifespan=lifespan, title="DivBase API", docs_url="/api/v1/docs")

app.include_router(users_router, prefix="/api/v1/users", tags=["users"])
app.include_router(projects_router, prefix="/api/v1/projects", tags=["projects"])
app.include_router(auth_router, prefix="/api/v1/auth", tags=["auth"])
app.include_router(admin_router, prefix="/api/v1/admin", tags=["admin"])

app.include_router(fr_auth_router, prefix="/auth", tags=["frontend", "auth"])
app.include_router(fr_admin_router, prefix="/admin", tags=["frontend", "admin"])
app.include_router(fr_core_router, prefix="", tags=["frontend"])
app.include_router(fr_profile_router, prefix="/profile", tags=["frontend", "profile"])


register_exception_handlers(app)


# TODO - move below routes into routes dir when ready.
@app.get("/health")
def health():
    return {"status": "ok"}


@app.get("/query/")
def get_jobs_by_user(user_name: str = "Default User"):
    """
    TODO: user_name would later be determined by the authentication system.
    """
    task_items = get_task_history()
    return task_items


@app.get("/query/{task_id}")
def get_task_by_id(task_id: str):
    task_items = get_task_history(task_id=task_id)
    return task_items


@app.post("/query/sample-metadata/")
def sample_metadata_query(tsv_filter: str, metadata_tsv_name: str, project: str):
    """
    Create a new bcftools query job for the specified project.
    """
    bucket_name = project
    task_kwargs = {
        "tsv_filter": tsv_filter,
        "metadata_tsv_name": metadata_tsv_name,
        "bucket_name": bucket_name,
    }

    results = sample_metadata_query_task.apply_async(kwargs=task_kwargs)
    result_dict = results.get(timeout=10)
    return result_dict


@app.post("/query/bcftools-pipe/")
def create_bcftools_jobs(
    tsv_filter: str,
    metadata_tsv_name: str,
    command: str,
    project: str,
    user_name: str = "Default User",
):
    """
    Create a new bcftools query job for the specified project.
    """
    bucket_name = project
    task_kwargs = {
        "tsv_filter": tsv_filter,
        "command": command,
        "metadata_tsv_name": metadata_tsv_name,
        "bucket_name": bucket_name,
        "user_name": user_name,
    }

    results = bcftools_pipe_task.apply_async(kwargs=task_kwargs)
    return results.id


@app.post("/dimensions/update/")
def update_vcf_dimensions_for_a_project(project: str, user_name: str = "Default User"):
    """
    Update the VCF dimensions files for the specified project
    """
    bucket_name = project
    task_kwargs = {
        "bucket_name": bucket_name,
        "user_name": user_name,
    }

    results = update_vcf_dimensions_task.apply_async(kwargs=task_kwargs)
    return results.id


def main():
    uvicorn.run("divbase_api.divbase_api:app", host="127.0.0.1", port=8000, reload=True)


if __name__ == "__main__":
    main()
