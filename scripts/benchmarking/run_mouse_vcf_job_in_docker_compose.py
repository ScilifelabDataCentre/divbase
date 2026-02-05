"""
A script that runs commands with divbase CLI and auxiliary scripts to perform the work-flow required to add a new VCF file to divbase (+ mock sample metadata)
and to run a bcftools query on the files. It downloads a large public VCF file and submits a divbase bcftools-query job.
The query (=celery task) is expected to take about 10 minutes to complete, depending on the machine.

NOTE! this script is intended to be run with the Docker Compose stack for local development. It has not been tested against a remote deployment of DivBase.

In addition to the query, the pre-processing will take a few minutes too
(downloading the file with curl, generating mock-metadata, check if there is a bucket called 'benchmarking',
uploading files to the local bucket if not already there).

In all, this little work-flow should be representative of all the steps needed to query a VCF file in the version of
Divbase at the time of the commit.

Usage:

    Ensure that divbase docker compose stack is running, e.g. with
    docker compose -f docker/divbase_compose.yaml up --build -d
    then initialize the local environment with:

    python scripts/benchmarking/local_dev_setup.py

    and then:

    python scripts/benchmarking/run_mouse_vcf_job.py
"""

import contextlib
import os
import re
import shlex
import subprocess
import sys
import time

import boto3
import httpx
from _benchmarking_shared_utils import LOCAL_ENV

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import local_dev_setup

local_dev_setup.PROJECTS = [
    {
        "name": "benchmarking-mouse",
        "description": "Benchmarking project",
        "bucket_name": "divbase-local-benchmarking-mouse",  # to comply with s3_service_account_policy.json
        "storage_quota_bytes": 10737418240,
        "files": ["README.md"],
    }
]
# they way this is implemented in local_dev_setup.py, files cannot be empty. So just upload README.md for now


def ensure_minio_bucket_exists(bucket_name: str) -> None:
    """Ensure the Minio bucket exists and has versioning enabled."""
    print(f"Ensuring Minio bucket '{bucket_name}' exists...")

    s3_client = boto3.client(
        "s3",
        endpoint_url=local_dev_setup.MINIO_URL,
        aws_access_key_id=local_dev_setup.MINIO_FAKE_ACCESS_KEY,
        aws_secret_access_key=local_dev_setup.MINIO_FAKE_SECRET_KEY,
    )

    # Try to create bucket, suppress exception (=exit) if it already exists
    with contextlib.suppress(s3_client.exceptions.BucketAlreadyOwnedByYou, s3_client.exceptions.BucketAlreadyExists):
        s3_client.create_bucket(Bucket=bucket_name)

    s3_client.put_bucket_versioning(
        Bucket=bucket_name,
        VersioningConfiguration={"Status": "Enabled"},
    )
    print(f"Bucket '{bucket_name}' is ready.")


def ensure_required_files_in_bucket(project_name: str, filename: str, mock_metadata: str, url: str) -> None:
    """Check if the required files are in the bucket. Otherwise: download the vcf file from the URL from the public repository;
    generate the mock metadata; upload the files to the bucket so that the worker container can access them from there during the query."""

    print("\nChecking if files are already in bucket...")
    cmd = f"divbase-cli files list --project {project_name}"
    env = dict(os.environ)
    env["DIVBASE_ENV"] = "local"

    result = subprocess.run(shlex.split(cmd), env=env, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        print(f"ERROR running 'divbase-cli files list': {result.stderr}")
        print(f"stdout: {result.stdout}")
        raise subprocess.CalledProcessError(result.returncode, cmd, output=result.stdout, stderr=result.stderr)
    output = result.stdout
    filenames_to_check = [filename, mock_metadata]

    for filename_to_check in filenames_to_check:
        if filename_to_check in output:
            print(f"{filename_to_check} is already present in the bucket.")
        else:
            print(
                f"{filename_to_check} is NOT present in the bucket. Checking to see if it exists in the current working directory..."
            )

            if not os.path.isfile(filename_to_check):
                print(f"{filename_to_check} not found in the current working directory.")

                if filename_to_check == filename:
                    cmd = f"curl -o {filename} {url}"
                    subprocess.run(shlex.split(cmd), check=True)
                    print("VCF file downloaded. Uploading to bucket...")
                elif filename_to_check == mock_metadata:
                    cmd = f"python scripts/generate_mock_sample_metadata.py --vcf {filename} --output {mock_metadata}"
                    subprocess.run(shlex.split(cmd), check=True)
                    print("Mock metadata generated. Uploading to bucket...")
            else:
                print(f"{filename} already exists. Uploading to bucket...")

            cmd_upload = shlex.split(f"divbase-cli files upload --project {project_name} {filename_to_check} ")
            subprocess.run(cmd_upload, check=True, env=LOCAL_ENV)


def wait_for_task_completion(
    job_id: str, project_name: str, admin_token: str, timeout: int = 600, poll_interval: int = 5
) -> bool:
    """
    Poll the task history API until the task completes successfully or fails.
    Uses direct API calls to get task status.
    Returns True if task completed successfully, raises exception otherwise.
    """
    print(f"Waiting for task {job_id} to complete...")
    start_time = time.time()

    while time.time() - start_time < timeout:
        try:
            response = httpx.get(
                f"{local_dev_setup.BASE_URL}/v1/task-history/projects/{project_name}",
                headers={"Authorization": f"Bearer {admin_token}"},
                params={"limit": 50},
            )
            response.raise_for_status()
            tasks = response.json()

            for task in tasks:
                if str(task.get("id")) == str(job_id):
                    status = task.get("status")
                    if status == "SUCCESS":
                        print(f"Task {job_id} completed successfully!")
                        return True
                    elif status == "FAILURE":
                        print(f"Task {job_id} failed!")
                        print(f"Result: {task.get('result')}")
                        raise RuntimeError(f"Task {job_id} failed")
                    else:
                        print(f"Task {job_id} status: {status} (elapsed: {int(time.time() - start_time)}s)")
                        break
            else:
                print(f"Task {job_id} not yet visible in history... (elapsed: {int(time.time() - start_time)}s)")

        except Exception as e:
            print(f"Error checking task status: {e}")

        time.sleep(poll_interval)

    raise TimeoutError(f"Task {job_id} did not complete within {timeout} seconds")


def ensure_project_exists_and_assign_manager(project_name: str, admin_token: str) -> None:
    """
    Wrapper around local_dev_setup.create_projects() that allows the script to be run multiple times.
    First time, the project will not exist and it will be created. Next time, it will see that the project
    exists and skip the creation step.
    """
    try:
        for project in local_dev_setup.PROJECTS:
            project_data = {
                "name": project["name"],
                "description": project["description"],
                "bucket_name": project["bucket_name"],
                "storage_quota_bytes": project["storage_quota_bytes"],
            }
            response = local_dev_setup.make_authenticated_request(
                "POST", f"{local_dev_setup.BASE_URL}/v1/admin/projects", admin_token, json=project_data
            )

            project = response.json()
            project_id = project["id"]
            print(f"Created project: {project_id}")

        local_dev_setup.make_authenticated_request(
            "POST",
            f"{local_dev_setup.BASE_URL}/v1/admin/projects/{project_id}/members/1",
            admin_token,
            params={"role": "manage"},
        )
        print(f"Assigned first admin user as manager to {project_name}")

    except Exception as e:
        if "already in use" in str(e):
            print(f"Project '{project_name}' already exists in the API. Skipping creation and role assignment.")
        else:
            raise


def main():
    filename = "mgp.v3.snps.rsIDdbSNPv137.vcf.gz"
    url = "ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ022/ERZ022025/mgp.v3.snps.rsIDdbSNPv137.vcf.gz"
    mock_metadata = "mock_metadata_mgpv3snps.tsv"
    project = local_dev_setup.PROJECTS[0]
    project_name = project["name"]
    bucket_name = project["bucket_name"]

    ensure_minio_bucket_exists(bucket_name)

    admin_token = local_dev_setup.get_admin_access_token()
    ensure_project_exists_and_assign_manager(project_name=project_name, admin_token=admin_token)

    local_dev_setup.create_local_config()

    # Login non-interactively by providing 'y' to the prompt, using credentials from local_dev_setup
    admin_email = local_dev_setup.ADMIN_CREDENTIALS["username"]
    admin_password = local_dev_setup.ADMIN_CREDENTIALS["password"]
    cmd = shlex.split(f"divbase-cli auth login {admin_email} --password {admin_password}")
    subprocess.run(cmd, input="y\n", text=True, check=True, env=LOCAL_ENV)

    ensure_required_files_in_bucket(project_name=project_name, filename=filename, mock_metadata=mock_metadata, url=url)

    # Run dimensions update and wait for completion
    print("\nUpdating VCF dimensions...")
    cmd_dimensions = f"divbase-cli dimensions update --project {project_name}"
    result = subprocess.run(shlex.split(cmd_dimensions), env=LOCAL_ENV, capture_output=True, text=True, check=True)

    # Extract job_id from the output
    # Expected format: "Job submitted successfully with task id: 102"
    match = re.search(r"task id:\s*(\d+)", result.stdout)
    if not match:
        print(f"Could not extract job_id from dimensions update command output:\n{result.stdout}")
        raise ValueError("Failed to extract job_id from dimensions update command")

    dimensions_job_id = match.group(1)
    print(f"Dimensions update job ID: {dimensions_job_id}")

    # Wait for dimensions update to complete
    wait_for_task_completion(dimensions_job_id, project_name=project_name, admin_token=admin_token, timeout=600)

    print("\nSubmitting query job to task queue...")
    cmd_query = f"divbase-cli query bcftools-pipe --tsv-filter 'Area:North,East' --command 'view -s SAMPLES; view -r 1:15000000-25000000' --metadata-tsv-name {mock_metadata} --project {project_name}"
    subprocess.run(shlex.split(cmd_query), check=True, env=LOCAL_ENV)

    print("\nThe query has been submitted. Check the divbase task history for updates.")


if __name__ == "__main__":
    main()
