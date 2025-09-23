"""
A helper script derived from local_dev_setup.py that is intended for use with a local k3d cluster.
Kept separate from local_dev_setup.py which uses the Docker Compose setup.


### Usage:
Once all pods are ready, you can run this script with:
uv run scripts/local_k3d_dev_setup.py # (or python scripts/local_k3d_dev_setup.py)
"""

import contextlib
import os
import shlex
import subprocess
from pathlib import Path

import boto3

PROJECTS = {
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


MINIO_URL = "http://localhost:30900"
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

    for project in PROJECTS:
        with contextlib.suppress(s3_client.exceptions.BucketAlreadyOwnedByYou):
            s3_client.create_bucket(Bucket=project)

        s3_client.put_bucket_versioning(
            Bucket=project,
            VersioningConfiguration={"Status": "Enabled"},
        )


def create_local_config():
    """Create a local user DivBase config file and add projects to it."""
    command = shlex.split("divbase-cli config create")
    result = subprocess.run(command, check=False, env=LOCAL_ENV, stderr=subprocess.PIPE)
    if result.returncode != 0 and "FileExistsError" not in result.stderr.decode():
        raise subprocess.CalledProcessError(result.returncode, command, output=result.stdout, stderr=result.stderr)

    for project_name in PROJECTS:
        command = shlex.split(
            f"divbase-cli config add-project {project_name} --divbase-url http://localhost:30800 --s3-url http://localhost:30900"
        )
        subprocess.run(command, check=True, env=LOCAL_ENV)

    default_project = list(PROJECTS.keys())[0]
    command = shlex.split(f"divbase-cli config set-default {default_project}")
    subprocess.run(command, check=True, env=LOCAL_ENV)


def upload_files_to_buckets():
    """Upload files to each created project's storage bucket."""
    for project_name, files in PROJECTS.items():
        file_list = [str(FIXTURES_DIR / file) for file in files]
        files_to_upload = " ".join(file_list)
        command = shlex.split(f"divbase-cli files upload --project {project_name} {files_to_upload}")
        subprocess.run(command, check=True, env=LOCAL_ENV)


def add_bucket_versioning_file():
    """Add a bucket versioning file to each project's storage bucket."""
    for project in PROJECTS:
        command = shlex.split(
            f"divbase-cli version add v0.1.0 --description 'add initial data sets' --project {project}"
        )
        subprocess.run(command, check=True, env=LOCAL_ENV)


if __name__ == "__main__":
    setup_minio_buckets()
    create_local_config()
    upload_files_to_buckets()
    add_bucket_versioning_file()
