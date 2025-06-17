"""
Helper module to set up and tear down a Minio instance used for testing purposes.
"""

import shlex
import subprocess
import time
from pathlib import Path

import boto3

URL = "http://localhost:9000"
MINIO_FAKE_ACCESS_KEY = "minioadmin"
MINIO_FAKE_SECRET_KEY = "badpassword"

DOCKER_COMPOSE_FILE = "docker/minio/docker-compose.yaml"
COMPOSE_COMMAND = shlex.split(f"docker-compose -f {DOCKER_COMPOSE_FILE} up -d minio")
HEALTH_CHECK_COMMAND = shlex.split(f"curl -I {URL}/minio/health/live")
STOP_COMMAND = shlex.split(f"docker-compose -f {DOCKER_COMPOSE_FILE} down")


BUCKETS = ["bucket1", "bucket2", "empty-bucket"]
FIXTURES_DIR = Path(__file__).parent.parent / "tests" / "fixtures"
UPLOADED_FILES = ["file1.txt", "file2.txt", "file3.txt"]


def start_minio() -> None:
    """Start Minio container using Docker compose, wait till available."""
    subprocess.run(COMPOSE_COMMAND, check=True)
    print("Waiting for Minio to start...")
    time.sleep(2)

    max_retries = 10
    for _ in range(max_retries):
        health_check = subprocess.run(HEALTH_CHECK_COMMAND, check=False, capture_output=True)
        if health_check.returncode == 0:
            return
        max_retries -= 1
        time.sleep(2)
        print("Still waiting for Minio to start...")
        continue

    raise RuntimeError("Minio failed to start properly within its allocated time.")


def stop_minio() -> None:
    subprocess.run(STOP_COMMAND, check=True)
    time.sleep(2)


def setup_minio_data() -> None:
    """Create test buckets and add some files to them in Minio."""
    s3_client = boto3.client(
        "s3",
        endpoint_url=URL,
        aws_access_key_id=MINIO_FAKE_ACCESS_KEY,
        aws_secret_access_key=MINIO_FAKE_SECRET_KEY,
    )

    for bucket in BUCKETS:
        s3_client.create_bucket(Bucket=bucket)

    for file in UPLOADED_FILES:
        s3_client.upload_file(Filename=str(FIXTURES_DIR / file), Bucket="bucket1", Key=file)
    s3_client.upload_file(Filename=str(FIXTURES_DIR / "file1.txt"), Bucket="bucket2", Key="file1.txt")
