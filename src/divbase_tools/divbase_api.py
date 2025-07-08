"""
The API server for DivBase.
"""

from pathlib import Path

import uvicorn
from fastapi import FastAPI

from divbase_tools.task_history import get_task_history
from divbase_tools.tasks import bcftools_full_query_task

TSV_FILE = Path("./sample_metadata.tsv")

app = FastAPI()


@app.get("/")
async def root():
    return {"message": "DivBase API is running!"}


@app.get("/health")
def health():
    return {"status": "ok"}


@app.post("/jobs/")
def create_job(tsv_filter: str, command: str, bucket_name: str, user_name: str = "Default User"):
    """
    Create a new query job in the specified bucket.

    TODO: user_name would later be determined by the authentication system.
    """
    task_kwargs = {
        "tsv_filter": tsv_filter,
        "command": command,
        "bucket_name": bucket_name,
        "user_name": user_name,
    }

    results = bcftools_full_query_task.apply_async(kwargs=task_kwargs)
    return results.id


@app.get("/jobs/")
def get_jobs_by_user(user_name: str = "Default User"):
    """
    TODO: user_name would later be determined by the authentication system.
    """
    task_items = get_task_history()
    return task_items


def main():
    uvicorn.run("divbase_tools.divbase_api:app", host="127.0.0.1", port=8000, reload=True)


if __name__ == "__main__":
    main()
