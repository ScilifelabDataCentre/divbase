"""
A helper script to set up some buckets and add them to user config for DivBase local development.
Can be used with both the local docker compose setup (docker/divbase_compose.yaml) and the
local k3d setup (stored in a separate private repository).

### This script will:
1. Create some buckets in Minio and populate them with test data from the fixtures directory.
2. Create a user config file for you and add these buckets to it.

### Requirements

#### Docker compose:
Have the local docker compose stack running and healthy:

docker compose -f docker/divbase_compose.yaml up -d

# to check if all healthy
docker compose -f docker/divbase_compose.yaml ps

# You could consider this one liner to first clean (DELETE ALL VOLUMES) and start fresh:
docker compose -f docker/divbase_compose.yaml down -v && docker compose -f docker/divbase_compose.yaml up --build -d

#### k3d:
Have the local k3d stack running and all pods ready as described in the private repository.

# to check if all healthy;
kubectl get pods

# the add the --k3d flag to the command below.

### Usage:
Once all containers are healthy, you can run this script with:
uv run scripts/local_dev_setup.py # (or python scripts/local_dev_setup.py)

Alternatively, add the --k3d flag if using the k3d setup, e.g.
python scripts/local_dev_setup.py --k3d

"""

import argparse
import contextlib
import os
import shlex
import subprocess
from pathlib import Path

import boto3

COMPOSE_PROJECTS = {
    "local-project-1": [
        "HOM_20ind_17SNPs_first_10_samples.vcf.gz",
        "HOM_20ind_17SNPs_last_10_samples.vcf.gz",
        "sample_metadata.tsv",
    ],
    "local-project-2": ["sample_metadata.tsv"],
    "local-project-3": [
        "HOM_20ind_17SNPs.1.vcf.gz",
        "HOM_20ind_17SNPs.13.vcf.gz",
        "HOM_20ind_17SNPs.18.vcf.gz",
        "HOM_20ind_17SNPs.20.vcf.gz",
        "HOM_20ind_17SNPs.21.vcf.gz",
        "HOM_20ind_17SNPs.22.vcf.gz",
        "HOM_20ind_17SNPs.24.vcf.gz",
        "HOM_20ind_17SNPs.4.vcf.gz",
        "HOM_20ind_17SNPs.5.vcf.gz",
        "HOM_20ind_17SNPs.6.vcf.gz",
        "HOM_20ind_17SNPs.7.vcf.gz",
        "HOM_20ind_17SNPs.8.vcf.gz",
        "sample_metadata_HOM_chr_split_version.tsv",
    ],
    "local-project-4": [
        "HOM_20ind_17SNPs.1.vcf.gz",
        "HOM_20ind_17SNPs.4.vcf.gz",
        "HOM_20ind_17SNPs.21.vcf.gz",
        "HOM_20ind_17SNPs.8_edit_new_sample_names.vcf.gz",
        "HOM_20ind_17SNPs.13_edit_new_sample_names.vcf.gz",
        "HOM_20ind_17SNPs.18_edit_new_sample_names.vcf.gz",
        "HOM_20ind_17SNPs_changed_sample_names.vcf.gz",
        "sample_metadata_HOM_files_that_need_mixed_bcftools_concat_and_merge.tsv",
    ],
}

K3D_PROJECTS = {
    "k3d-test-project-1": [
        "HOM_20ind_17SNPs_first_10_samples.vcf.gz",
        "HOM_20ind_17SNPs_last_10_samples.vcf.gz",
        "sample_metadata.tsv",
    ],
    "k3d-test-project-2": ["sample_metadata.tsv"],
    "k3d-test-project-3": [
        "HOM_20ind_17SNPs.1.vcf.gz",
        "HOM_20ind_17SNPs.13.vcf.gz",
        "HOM_20ind_17SNPs.18.vcf.gz",
        "HOM_20ind_17SNPs.20.vcf.gz",
        "HOM_20ind_17SNPs.21.vcf.gz",
        "HOM_20ind_17SNPs.22.vcf.gz",
        "HOM_20ind_17SNPs.24.vcf.gz",
        "HOM_20ind_17SNPs.4.vcf.gz",
        "HOM_20ind_17SNPs.5.vcf.gz",
        "HOM_20ind_17SNPs.6.vcf.gz",
        "HOM_20ind_17SNPs.7.vcf.gz",
        "HOM_20ind_17SNPs.8.vcf.gz",
        "sample_metadata_HOM_chr_split_version.tsv",
    ],
    "k3d-test-project-4": [
        "HOM_20ind_17SNPs.1.vcf.gz",
        "HOM_20ind_17SNPs.4.vcf.gz",
        "HOM_20ind_17SNPs.21.vcf.gz",
        "HOM_20ind_17SNPs.8_edit_new_sample_names.vcf.gz",
        "HOM_20ind_17SNPs.13_edit_new_sample_names.vcf.gz",
        "HOM_20ind_17SNPs.18_edit_new_sample_names.vcf.gz",
        "HOM_20ind_17SNPs_changed_sample_names.vcf.gz",
        "sample_metadata_HOM_files_that_need_mixed_bcftools_concat_and_merge.tsv",
    ],
}

MINIO_FAKE_ACCESS_KEY = "minioadmin"
MINIO_FAKE_SECRET_KEY = "badpassword"

FIXTURES_DIR = Path(__file__).parent.parent / "tests" / "fixtures"

LOCAL_ENV = os.environ.copy()
LOCAL_ENV["DIVBASE_ENV"] = "local"


def setup_minio_buckets(projects, minio_url, access_key, secret_key):
    """Create buckets and enable bucket versioning"""
    s3_client = boto3.client(
        "s3",
        endpoint_url=minio_url,
        aws_access_key_id=access_key,
        aws_secret_access_key=secret_key,
    )

    for project in projects:
        with contextlib.suppress(s3_client.exceptions.BucketAlreadyOwnedByYou):
            s3_client.create_bucket(Bucket=project)

        s3_client.put_bucket_versioning(
            Bucket=project,
            VersioningConfiguration={"Status": "Enabled"},
        )


def create_local_config(projects, local_env, api_url, s3_url):
    """Create a local user DivBase config file and add projects to it."""
    command = shlex.split("divbase-cli config create")
    result = subprocess.run(command, check=False, env=local_env, stderr=subprocess.PIPE)
    if result.returncode != 0 and "FileExistsError" not in result.stderr.decode():
        raise subprocess.CalledProcessError(result.returncode, command, output=result.stdout, stderr=result.stderr)

    for project_name in projects:
        command = shlex.split(
            f"divbase-cli config add-project {project_name} --divbase-url {api_url} --s3-url {s3_url}"
        )
        subprocess.run(command, check=True, env=local_env)

    default_project = list(projects.keys())[0]
    command = shlex.split(f"divbase-cli config set-default {default_project}")
    subprocess.run(command, check=True, env=local_env)


def upload_files_to_buckets(projects, local_env):
    """Upload files to each created project's storage bucket."""
    for project_name, files in projects.items():
        file_list = [str(FIXTURES_DIR / file) for file in files]
        files_to_upload = " ".join(file_list)
        command = shlex.split(f"divbase-cli files upload --project {project_name} {files_to_upload}")
        subprocess.run(command, check=True, env=local_env)


def add_bucket_versioning_file(projects, local_env):
    """Add a bucket versioning file to each project's storage bucket."""
    for project in projects:
        command = shlex.split(
            f"divbase-cli version add v0.1.0 --description 'add initial data sets' --project {project}"
        )
        subprocess.run(command, check=True, env=local_env)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="DivBase local/k3d dev setup")
    parser.add_argument("--k3d", action="store_true", help="Use k3d cluster settings")
    args = parser.parse_args()

    MINIO_URL = "http://localhost:9000"
    API_URL = "http://localhost:8000"
    S3_URL = "http://localhost:9000"

    if args.k3d:
        PROJECTS = K3D_PROJECTS
        ## for when different ports are needed for k3d
        # MINIO_URL = "http://localhost:30900"
        # API_URL = "http://localhost:30800"
        # S3_URL = "http://localhost:30900"
        ENV = "local-k3d"
    else:
        PROJECTS = COMPOSE_PROJECTS
        ENV = "local"

    setup_minio_buckets(PROJECTS, MINIO_URL, "minioadmin", "badpassword")
    create_local_config(PROJECTS, LOCAL_ENV, API_URL, S3_URL)
    upload_files_to_buckets(PROJECTS, LOCAL_ENV)
    add_bucket_versioning_file(PROJECTS, LOCAL_ENV)
