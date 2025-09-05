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
    "split-scaffold-project": [
        "HOM_20ind_17SNPs.1.vcf.gz",
        "HOM_20ind_17SNPs.4.vcf.gz",
        "HOM_20ind_17SNPs.5.vcf.gz",
        "HOM_20ind_17SNPs.6.vcf.gz",
        "HOM_20ind_17SNPs.7.vcf.gz",
        "HOM_20ind_17SNPs.8.vcf.gz",
        "HOM_20ind_17SNPs.13.vcf.gz",
        "HOM_20ind_17SNPs.18.vcf.gz",
        "HOM_20ind_17SNPs.20.vcf.gz",
        "HOM_20ind_17SNPs.21.vcf.gz",
        "HOM_20ind_17SNPs.22.vcf.gz",
        "HOM_20ind_17SNPs.24.vcf.gz",
        "sample_metadata_HOM_chr_split_version.tsv",
    ],
    "cleaned-project": [],  # this project's bucket is always cleaned before a test
    "empty-project": [],
    "mixed-concat-merge-project": [
        "HOM_20ind_17SNPs.1.vcf.gz",
        "HOM_20ind_17SNPs.4.vcf.gz",
        "HOM_20ind_17SNPs.21.vcf.gz",
        "HOM_20ind_17SNPs_changed_sample_names.vcf.gz",
        "sample_metadata_HOM_files_that_need_mixed_bcftools_concat_and_merge.tsv",
    ],
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
