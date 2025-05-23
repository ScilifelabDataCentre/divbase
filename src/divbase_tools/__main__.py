"""
Manage versions of the entire bucket

This script is interacted with via argparse.

Versioning specified by updating a file in the bucket called ".bucket_versions.yaml".

Versioning done by setting the version according to the current timestamp.
The version file is always kept in the bucket, not stored on disk.

TODO:
- Check how robust timestamping approach is.
- Add support to delete versions?
- Add support to specify a timestamp for a version?
"""

from pathlib import Path

import typer

from divbase_tools.s3_client import S3FileManager, get_credentials
from divbase_tools.services import (
    add_version_command,
    create_version_object_command,
    download_files_command,
    list_versions_command,
)

MINIO_URL = "api.divbase-testground.scilifelab-2-dev.sys.kth.se"

app = typer.Typer(
    help="""This tool lets you interact with your DivBase project's bucket and do:
        - Upload and download files to/from the bucket.
        - Managed the Version the state of the entire bucket.
    """,
    no_args_is_help=True,
    context_settings={"help_option_names": ["-h", "--help"]},
)

file_app = typer.Typer()
version_app = typer.Typer()


@app.callback()
def main(
    ctx: typer.Context,
):
    """This creates the s3 file manager which is passed into all commands."""
    access_key, secret_key = get_credentials()

    ctx.obj = S3FileManager(
        url=MINIO_URL,
        access_key=access_key,
        secret_key=secret_key,
    )


@version_app.command("create")
def create(
    ctx: typer.Context, bucket_name: str = typer.Option(..., help="Name of the storage bucket for the project.")
):
    """Create the versioning metadata file to store the bucket in."""
    create_version_object_command(bucket_name, s3_file_manager=ctx.obj)


@version_app.command("add")
def add(
    ctx: typer.Context,
    name: str = typer.Option(..., help="Name of the version (e.g., semantic version)."),
    description: str = typer.Option(None, help="Optional description of the version."),
    bucket_name: str = typer.Option(..., help="Name of the storage bucket for the project."),
):
    """Add a new bucket version."""
    add_version_command(bucket_name, name, description, s3_file_manager=ctx.obj)


@version_app.command("list")
def list(
    ctx: typer.Context,
    bucket_name: str = typer.Option(..., help="Name of the storage bucket for the project."),
):
    """List all bucket versions."""
    list_versions_command(bucket_name, s3_file_manager=ctx.obj)


@version_app.command("delete")
def version_delete():
    pass


@file_app.command("download")
def download_files(
    ctx: typer.Context,
    bucket_name: str = typer.Option(..., help="Name of the storage bucket for the project."),
    files: str = typer.Option(..., help="Comma separated list of files to download from the bucket."),
    download_dir: Path = typer.Option(..., help="Directory to download the files to."),
    bucket_version: str = typer.Option(None, help="Version of the bucket at which to download the files."),
):
    """Download files from the bucket."""
    if not download_dir:
        download_dir = Path(".")

    download_files_command(
        bucket_name=bucket_name,
        files=files,
        download_dir=download_dir,
        bucket_version=bucket_version,
        s3_file_manager=ctx.obj,
    )


@file_app.command("upload")
def upload_files():
    pass


app.add_typer(version_app, name="version")
app.add_typer(file_app, name="file")


if __name__ == "__main__":
    app()
