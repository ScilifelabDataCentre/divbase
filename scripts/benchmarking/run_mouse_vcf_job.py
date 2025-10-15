"""
A script that runs commands with divbase CLI and auxiliary scripts to perform the work-flow required to add a new VCF file to divbase (+ mock sample metadata)
and to run a bcftools query on the files. It downloads a large public VCF file and submits a divbase bcftools-query job.
The query (=celery task) is expected to take about 10 minutes to complete, depending on the machine.

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

import os
import shlex
import subprocess
import sys

import boto3
import yaml
from _benchmarking_shared_utils import LOCAL_ENV

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import local_dev_setup

local_dev_setup.PROJECTS = [
    {
        "name": "benchmarking-mouse",
        "description": "Benchmarking project",
        "bucket_name": "benchmarking-mouse",
        "storage_quota_bytes": 10737418240,
        "files": ["README.md"],
    }
]
# they way this is implemented in local_dev_setup.py, files cannot be empty. So just upload README.md for now


def ensure_required_files_in_bucket(project_name: str, filename: str, mock_metadata: str, url: str) -> None:
    """Check if the required files are in the bucket. Otherwise: download the vcf file from the URL from the public repository;
    generate the mock metadata; upload the files to the bucket so that the worker container can access them from there during the query."""

    print("\nChecking if files are already in bucket...")
    cmd = f"divbase-cli files list --project {project_name}"
    env = dict(os.environ)
    env["DIVBASE_ENV"] = "local"

    result = subprocess.run(shlex.split(cmd), env=env, capture_output=True, text=True, check=True)
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


def upload_mock_vcf_dimensions_file(project_name: str, filename: str) -> None:
    """
    Upload a pre-calculated dimensions file to the bucket to save time (and avoid timing issues since the query task requires the file).
    Need to patch the fixture with the version number of the latest version of the VCF file in the bucket.
    """

    s3 = boto3.client(
        "s3",
        endpoint_url=local_dev_setup.MINIO_URL,
        aws_access_key_id=local_dev_setup.MINIO_FAKE_ACCESS_KEY,
        aws_secret_access_key=local_dev_setup.MINIO_FAKE_SECRET_KEY,
    )

    response = s3.list_object_versions(Bucket=project_name, Prefix=filename)
    latest_version = None
    for version in response.get("Versions", []):
        if version.get("IsLatest"):
            latest_version = version["VersionId"]
            break

    # Check if there is a .vcf_dimensions.yaml in the bucket already, and which VCF file version it refers to
    try:
        response = s3.get_object(Bucket=project_name, Key=".vcf_dimensions.yaml")
        s3_yaml_str = response["Body"].read().decode("utf-8")
    except s3.exceptions.NoSuchKey:
        s3_yaml_str = None

    if s3_yaml_str:
        dimensions_data = yaml.safe_load(s3_yaml_str)
        for entry in dimensions_data.get("dimensions", []):
            if (
                entry.get("filename") == "mgp.v3.snps.rsIDdbSNPv137.vcf.gz"
                and entry.get("file_version_ID_in_bucket") == latest_version
            ):
                print("Bucket already contains .vcf_dimensions.yaml with the latest version ID. Skipping upload.")
                return

    # If not latest VCF version  .vcf_dimensions.yaml in the bucket in update and upload a new version
    with open("tests/fixtures/vcf_dimensions_mgp.v3.snps.rsIDdbSNPv137.yaml", "r") as f:
        dimensions_data = yaml.safe_load(f)

    for entry in dimensions_data.get("dimensions", []):
        if entry.get("filename") == "mgp.v3.snps.rsIDdbSNPv137.vcf.gz":
            entry["file_version_ID_in_bucket"] = latest_version

    with open(".vcf_dimensions.yaml", "w") as f:
        yaml.safe_dump(dimensions_data, f)

    cmd_upload = shlex.split(f"divbase-cli files upload --project {project_name} .vcf_dimensions.yaml")
    subprocess.run(cmd_upload, check=True, env=LOCAL_ENV)
    os.remove(".vcf_dimensions.yaml")
    print("Uploaded new .vcf_dimensions.yaml with updated version ID.")


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

    admin_token = local_dev_setup.get_admin_access_token()
    ensure_project_exists_and_assign_manager(project_name=project_name, admin_token=admin_token)

    local_dev_setup.create_local_config()

    local_dev_setup.setup_minio_buckets()

    ensure_required_files_in_bucket(project_name=project_name, filename=filename, mock_metadata=mock_metadata, url=url)

    upload_mock_vcf_dimensions_file(project_name=project_name, filename=filename)

    print("\nSubmitting query job to task queue...")
    cmd_query = f"divbase-cli query bcftools-pipe --tsv-filter 'Area:North,East' --command 'view -s SAMPLES; view -r 1:15000000-25000000' --metadata-tsv-name {mock_metadata} --project {project_name}"
    subprocess.run(shlex.split(cmd_query), check=True, env=LOCAL_ENV)

    print("\nThe query has been submitted. Check the task status or the flower logs for updates.")


if __name__ == "__main__":
    main()
