"""
Manage the starlette-admin interface for the DivBase API.
"""

from fastapi import FastAPI
from sqlalchemy.ext.asyncio import AsyncEngine
from starlette_admin.contrib.sqla import Admin, ModelView

from divbase_api.models.projects import ProjectDB, ProjectMembershipDB
from divbase_api.models.users import UserDB


def register_admin_panel(app: FastAPI, engine: AsyncEngine) -> None:
    """
    Create and register an admin panel to the FastAPI app.
    """
    admin = Admin(engine, title="DivBase Admin")

    admin.add_view(ModelView(UserDB, icon="fas fa-user"))
    admin.add_view(ModelView(ProjectDB, icon="fas fa-project"))
    admin.add_view(ModelView(ProjectMembershipDB, icon="fas fa-link"))

    admin.mount_to(app)
