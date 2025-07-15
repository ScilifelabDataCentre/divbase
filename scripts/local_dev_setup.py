"""
A helper script to set up some buckets and add them to user config for DivBase local development.

### This script will:
1. Create some buckets in Minio and populate them with test data from the fixtures directory.
2. Create a user config file for you and add these buckets to it.

### Requirements

1. Have the local docker compose stack running and healthy:

docker compose -f docker/divbase_compose.yaml up -d

# to check if all healthy
docker compose -f docker/divbase_compose.yaml ps

# You could consider this one liner to first clean (DELETE ALL VOLUMES) and start fresh:
docker compose -f docker/divbase_compose.yaml down -v && docker compose -f docker/divbase_compose.yaml up --build -d

### Usage:
Once all containers are healthy, you can run this script with:
uv run scripts/local_dev_setup.py # (or python scripts/local_dev_setup.py)
"""

import contextlib
import os
import shlex
import subprocess
from pathlib import Path

import boto3

BUCKETS = {
    "local-bucket": [
        "HOM_20ind_17SNPs_first_10_samples.vcf.gz",
        "HOM_20ind_17SNPs_last_10_samples.vcf.gz",
        "sample_metadata.tsv",
    ],
    "local-bucket-2": ["sample_metadata.tsv"],
}

MINIO_URL = "http://localhost:9000"
MINIO_FAKE_ACCESS_KEY = "minioadmin"
MINIO_FAKE_SECRET_KEY = "badpassword"

FIXTURES_DIR = Path(__file__).parent.parent / "tests" / "fixtures"

LOCAL_ENV = os.environ.copy()
LOCAL_ENV["DIVBASE_ENV"] = "local"


def setup_minio_buckets() -> None:
    """Create buckets and enable bucket versioning"""
    s3_client = boto3.client(
        "s3",
        endpoint_url=MINIO_URL,
        aws_access_key_id=MINIO_FAKE_ACCESS_KEY,
        aws_secret_access_key=MINIO_FAKE_SECRET_KEY,
    )

    for bucket in BUCKETS:
        with contextlib.suppress(s3_client.exceptions.BucketAlreadyOwnedByYou):
            s3_client.create_bucket(Bucket=bucket)

        s3_client.put_bucket_versioning(
            Bucket=bucket,
            VersioningConfiguration={"Status": "Enabled"},
        )


def create_local_config():
    """Create a local user DivBase config file and add buckets to it."""
    command = shlex.split("divbase-cli config create")
    result = subprocess.run(command, check=False, env=LOCAL_ENV, stderr=subprocess.PIPE)
    if result.returncode != 0 and "FileExistsError" not in result.stderr.decode():
        raise subprocess.CalledProcessError(result.returncode, command, output=result.stdout, stderr=result.stderr)

    for bucket_name in BUCKETS:
        command = shlex.split(f"divbase-cli config add-bucket {bucket_name}")
        subprocess.run(command, check=True, env=LOCAL_ENV)

    default_bucket = list(BUCKETS.keys())[0]
    command = shlex.split(f"divbase-cli config set-default {default_bucket}")
    subprocess.run(command, check=True, env=LOCAL_ENV)


def upload_files_to_buckets():
    """Upload files to the created buckets."""
    for bucket_name, files in BUCKETS.items():
        file_list = [str(FIXTURES_DIR / file) for file in files]
        files_to_upload = " ".join(file_list)
        command = shlex.split(f"divbase-cli files upload --bucket-name {bucket_name} {files_to_upload}")
        subprocess.run(command, check=True, env=LOCAL_ENV)


def add_bucket_versioning_file():
    """Add a bucket versioning file to each bucket."""
    for bucket_name in BUCKETS:
        command = shlex.split(
            f"divbase-cli version add v0.1.0 --description 'add initial data sets'  --bucket-name {bucket_name}"
        )
        subprocess.run(command, check=True, env=LOCAL_ENV)


if __name__ == "__main__":
    setup_minio_buckets()
    create_local_config()
    upload_files_to_buckets()
    add_bucket_versioning_file()
