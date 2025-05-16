"""
This script is used to upload files to a MinIO bucket for a versioning experiment.

The plan is to version the entire bucket by using a yaml file called ".bucket_version.yaml" stored in the bucket.
The bucket level versioning is done by updating this file according to the current timestamp.

Note: Could not get Minio python client to work with user access key and secret key
(even though worked with: console, mc cli and boto3 python client)
and minio python client worked with admin credentials.

Therefore have decided to use boto3 python client for these experiments, issue can be revisited later.
"""

import os
from pathlib import Path

import boto3
import botocore
from dotenv import load_dotenv

MINIO_URL = "api.divbase-testground.scilifelab-2-dev.sys.kth.se"
BUCKET_NAME = "version-bucket-state"
SAMPLE_FILES = Path("experiments/version_entire_bucket/sample_files")
DOWNLOADS_DIR = Path("experiments/version_entire_bucket/downloads")


def list_objects(s3_client: "botocore.client.S3") -> None:
    response = s3_client.list_objects_v2(Bucket=BUCKET_NAME)
    if "Contents" in response:
        print(f"Objects in bucket '{BUCKET_NAME}':")
        for obj in response["Contents"]:
            print(obj["Key"])
    else:
        print(f"No objects found in bucket '{BUCKET_NAME}'.")


def download_file(s3_client: "botocore.client.S3", object_name: str, file_path: Path):
    s3_client.download_file(Bucket=BUCKET_NAME, Key=object_name, Filename=str(file_path))
    print(f"Object '{object_name}' downloaded to '{file_path}'.")


def upload_file(s3_client: "botocore.client.S3", object_name: str, file_path: Path):
    _ = s3_client.upload_file(Filename=str(file_path), Bucket=BUCKET_NAME, Key=object_name)
    print(f"File '{file_path}' uploaded as '{object_name}'.")


if __name__ == "__main__":
    load_dotenv()
    ACCESS_KEY = os.getenv("MINIO_SQUIRREL_USER_ACCESS_KEY")
    SECRET_KEY = os.getenv("MINIO_SQUIRREL_USER_SECRET_KEY")

    s3_client = boto3.client(
        "s3",
        endpoint_url=f"https://{MINIO_URL}",
        aws_access_key_id=ACCESS_KEY,
        aws_secret_access_key=SECRET_KEY,
    )

    list_objects(s3_client=s3_client)

    # download_file_path = DOWNLOADS_DIR / "test_pic.png"
    # download_file(s3_client=s3_client, object_name="test_pic.png", file_path=download_file_path)

    # Object name is the name of the file in the bucket,
    # So you can upload 2 different files but set same object name to create a new version of the file.
    # NOTE: The below section is commented out so we don't run this again, by accident.

    # ori_file1_path = SAMPLE_FILES / "file1_original.txt"
    # upload_file(s3_client=s3_client, object_name="file1.txt", file_path=ori_file1_path)
    # file2_path = SAMPLE_FILES / "file2.txt"
    # upload_file(s3_client=s3_client, object_name="file2.txt", file_path=file2_path)

    # A new version v0.1.0 of the bucket was added between the uploads above and below.

    # file1_new_path = SAMPLE_FILES / "file1_new.txt"
    # upload_file(s3_client=s3_client, object_name="file1.txt", file_path=file1_new_path)
    # file3_path = SAMPLE_FILES / "file3.txt"
    # upload_file(s3_client=s3_client, object_name="file3.txt", file_path=file3_path)
