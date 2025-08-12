"""
A script that downloads a large public VCF file and submits a divbase bcftools-query job.
The query (=celery task) is expected to take about 10 minutes to complete, depending on the machine.

In addition to the query, the pre-processing will take a few minutes too
(downloading the file with curl, generating mock-metadata, check if there is a bucket called 'benchmarking',
uploading files to the local bucket if not already there).

In all, this little work-flow should be representative of all the steps needed to query a VCF file in the version of
Divbase at the time of the commit.

Usage:
    python scripts/run_mouse_vcf_job.py
"""

import contextlib
import os
import shlex
import subprocess
from pathlib import Path

import boto3
import yaml
from local_dev_setup import LOCAL_ENV, MINIO_FAKE_ACCESS_KEY, MINIO_FAKE_SECRET_KEY, MINIO_URL


def ensure_project_exists(project_name: str):
    """Create a Minio bucket and DivBase project config entry if they do not already exist."""

    print(f"\nChecking if there is a bucket named {project_name} in the local MinIO container...")

    s3_client = boto3.client(
        "s3",
        endpoint_url=MINIO_URL,
        aws_access_key_id=MINIO_FAKE_ACCESS_KEY,
        aws_secret_access_key=MINIO_FAKE_SECRET_KEY,
    )

    existing_buckets = [b["Name"] for b in s3_client.list_buckets().get("Buckets", [])]
    if project_name in existing_buckets:
        print(f"Bucket '{project_name}' already exists.")
    else:
        print(f"Bucket '{project_name}' does not exist. Creating it...")
        with contextlib.suppress(s3_client.exceptions.BucketAlreadyOwnedByYou):
            s3_client.create_bucket(Bucket=project_name)
        print(f"Bucket '{project_name}' created.")

    s3_client.put_bucket_versioning(
        Bucket=project_name,
        VersioningConfiguration={"Status": "Enabled"},
    )


def ensure_project_in_config(project_name: str) -> None:
    """Ensure that the project is listed in the DivBase config."""

    config_path = Path.home() / ".config" / ".divbase_tools.yaml"

    print(f"\nChecking if {project_name} exists in config...")

    if project_in_config(project_name, config_path):
        print(f"Project '{project_name}' already exists in config.")
    else:
        print(f"Project '{project_name}' does not exist in config. Adding it...")
        command = shlex.split(f"divbase-cli config add-project {project_name}")
        subprocess.run(command, check=True, env=LOCAL_ENV)


def project_in_config(project_name: str, config_path: Path) -> bool:
    if not config_path.exists():
        return False
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    projects = config.get("projects", [])
    return any(p.get("name") == project_name for p in projects)


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


def main():
    filename = "mgp.v3.snps.rsIDdbSNPv137.vcf.gz"
    url = "ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ022/ERZ022025/mgp.v3.snps.rsIDdbSNPv137.vcf.gz"
    mock_metadata = "mock_metadata_mgpv3snps.tsv"
    project_name = "benchmarking"

    ensure_project_exists(project_name=project_name)

    ensure_project_in_config(project_name=project_name)

    ensure_required_files_in_bucket(project_name=project_name, filename=filename, mock_metadata=mock_metadata, url=url)

    print("\nSubmitting query job to task queue...")
    cmd_query = f"divbase-cli query bcftools-pipe --tsv-filter 'Area:North,East' --command 'view -s SAMPLES; view -r 1:15000000-25000000' --metadata-tsv-name {mock_metadata} --project benchmarking"
    subprocess.run(shlex.split(cmd_query), check=True, env=LOCAL_ENV)

    print("\nThe query has been submitted. Check the task status or the flower logs for updates.")

    print(
        "\nCalculating dimensions of the VCF file... (Run as the final step of the script since it can take some time for large VCF files)\n"
    )
    cmd_calculate_dimensions = f"python scripts/calculate_dimensions_of_vcf.py --vcf {filename}"
    subprocess.run(shlex.split(cmd_calculate_dimensions), check=True)


if __name__ == "__main__":
    main()
