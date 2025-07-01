"""
The API server for DivBase.
"""

from pathlib import Path

import uvicorn
from fastapi import FastAPI

from divbase_tools.queries import SidecarQueryManager
from divbase_tools.tasks import bcftools_pipe_task

TSV_FILE = Path("./sample_metadata.tsv")

app = FastAPI()


@app.get("/")
async def root():
    return {"message": "DivBase API is running!"}


@app.post("/jobs/")
def create_job(tsv_filter: str, command: str, bucket_name: str, user_name: str = "Default User"):
    """
    Create a new query job in the specified bucket.

    user_name would later be determined by the authentication system.
    """
    sidecar_manager = SidecarQueryManager(file=TSV_FILE).run_query(filter_string=tsv_filter)

    unique_sampleIDs = sidecar_manager.get_unique_values("Sample_ID")
    unique_filenames = sidecar_manager.get_unique_values("Filename")
    sample_and_filename_subset = sidecar_manager.query_result[["Sample_ID", "Filename"]]
    serialized_samples = sample_and_filename_subset.to_dict(orient="records")

    bcftools_inputs = {
        "sample_and_filename_subset": serialized_samples,
        "sampleIDs": unique_sampleIDs,
        "filenames": unique_filenames,
    }
    # TODO - assuming all files needed already available locally.

    task_kwargs = {
        "command": command,
        "bcftools_inputs": bcftools_inputs,
        "submitter": user_name,
    }

    result = bcftools_pipe_task.apply_async(kwargs=task_kwargs)
    return result.id


@app.get("/jobs")
def get_jobs_by_user():
    """
    Retrieve a list of jobs submitted by the user.
    """
    return {"message": "List of jobs"}


def main():
    uvicorn.run("divbase_tools.divbase_api:app", host="127.0.0.1", port=8000, reload=True)


if __name__ == "__main__":
    main()
