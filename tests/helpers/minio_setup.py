"""
Helper module to set up and tear down a Minio instance used for testing purposes.
"""

from pathlib import Path

import boto3

MINIO_URL = "http://localhost:9002"
MINIO_FAKE_ACCESS_KEY = "minioadmin"
MINIO_FAKE_SECRET_KEY = "badpassword"

FIXTURES_DIR = Path(__file__).parent.parent / "fixtures"

PROJECTS = {
    "project1": ["file1.txt", "file2.txt"],
    "project2": ["file1.txt"],
    "query-project": [
        "HOM_20ind_17SNPs_first_10_samples.vcf.gz",
        "HOM_20ind_17SNPs_last_10_samples.vcf.gz",
        "sample_metadata.tsv",
    ],
    "cleaned-project": [],  # this project's bucket is always cleaned before a test
    "empty-project": [],
}


def setup_minio_data() -> None:
    """Create test buckets and add some files to them in Minio."""
    s3_client = boto3.client(
        "s3",
        endpoint_url=MINIO_URL,
        aws_access_key_id=MINIO_FAKE_ACCESS_KEY,
        aws_secret_access_key=MINIO_FAKE_SECRET_KEY,
    )

    for bucket_name in PROJECTS:
        s3_client.create_bucket(Bucket=bucket_name)
        s3_client.put_bucket_versioning(
            Bucket=bucket_name,
            VersioningConfiguration={"Status": "Enabled"},
        )

    for bucket_name, files in PROJECTS.items():
        for file in files:
            s3_client.upload_file(Filename=str(FIXTURES_DIR / file), Bucket=bucket_name, Key=file)
