"""
Routes for project management.
"""

from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.ext.asyncio import AsyncSession

from divbase_api.crud.projects import create_project, get_all_projects, get_project_by_id
from divbase_api.db import get_db
from divbase_api.schemas.projects import ProjectCreate, ProjectResponse

projects_router = APIRouter()


@projects_router.get("/", response_model=list[ProjectResponse], status_code=status.HTTP_200_OK)
async def get_all_projects_endpoint(db: AsyncSession = Depends(get_db), limit: int = 1000):
    project = await get_all_projects(db, limit=limit)
    return project


@projects_router.post("/", response_model=ProjectResponse, status_code=status.HTTP_201_CREATED)
async def create_project_endpoint(proj_data: ProjectCreate, db: AsyncSession = Depends(get_db)):
    new_project = await create_project(db=db, proj_data=proj_data)
    return new_project


@projects_router.get("/id/{project_id}", response_model=ProjectResponse, status_code=status.HTTP_200_OK)
async def get_project_by_id_endpoint(project_id: int, db: AsyncSession = Depends(get_db)):
    project = await get_project_by_id(db=db, id=project_id)
    if not project:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Project not found")
    return project


# TODO
@projects_router.post("/{project_id}/members/{user_id}")
async def add_project_member():
    pass


# TODO
@projects_router.get("/{project_id}/members")
async def get_project_members():
    pass
