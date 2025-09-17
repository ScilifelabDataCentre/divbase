"""
The API server for DivBase.

TODO: user_name would later be determined by the authentication system.
"""

from pathlib import Path

import uvicorn
from fastapi import FastAPI

from divbase_api.get_task_history import get_task_history
from divbase_worker.tasks import bcftools_pipe_task, sample_metadata_query_task, update_vcf_dimensions_task

TSV_FILE = Path("./sample_metadata.tsv")

app = FastAPI()


@app.get("/")
async def root():
    return {"message": "DivBase API is running"}


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
    tsv_filter: str, metadata_tsv_name: str, command: str, project: str, user_name: str = "Default User"
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
