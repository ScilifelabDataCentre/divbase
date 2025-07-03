"""
Helper module to set up and tear down a Minio instance used for testing purposes.
"""

from pathlib import Path

import boto3

MINIO_URL = "http://localhost:9002"  # from overide in docker compose tests.yml file
MINIO_FAKE_ACCESS_KEY = "minioadmin"
MINIO_FAKE_SECRET_KEY = "badpassword"

FIXTURES_DIR = Path(__file__).parent.parent / "fixtures"

BUCKETS = {
    "bucket1": ["file1.txt", "file2.txt"],
    "bucket2": ["file1.txt"],
    "cleaned-bucket": [],  # this bucket is always cleaned before a test
    "empty-bucket": [],
}


def setup_minio_data() -> None:
    """Create test buckets and add some files to them in Minio."""
    s3_client = boto3.client(
        "s3",
        endpoint_url=MINIO_URL,
        aws_access_key_id=MINIO_FAKE_ACCESS_KEY,
        aws_secret_access_key=MINIO_FAKE_SECRET_KEY,
    )

    for bucket in BUCKETS:
        s3_client.create_bucket(Bucket=bucket)
        s3_client.put_bucket_versioning(
            Bucket=bucket,
            VersioningConfiguration={"Status": "Enabled"},
        )

    for bucket, files in BUCKETS.items():
        for file in files:
            s3_client.upload_file(Filename=str(FIXTURES_DIR / file), Bucket=bucket, Key=file)
