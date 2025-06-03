from pathlib import Path
from typing import List

import typer
from rich import print
from typing_extensions import Annotated

from divbase_tools.cli.user_config_cli import CONFIG_PATH_OPTION
from divbase_tools.cli.utils import resolve_bucket_name
from divbase_tools.cli.version_cli import BUCKET_NAME_OPTION
from divbase_tools.services import (
    download_files_command,
    list_files_command,
    upload_files_command,
)

file_app = typer.Typer(no_args_is_help=True, help="Download/upload files to/from the bucket.")


@file_app.command("list")
def list_files(
    bucket_name: str = BUCKET_NAME_OPTION,
    config_path: Path = CONFIG_PATH_OPTION,
):
    """list all files in the bucket."""
    bucket_name = resolve_bucket_name(bucket_name, config_path)
    files = list_files_command(
        bucket_name=bucket_name,
        config_path=config_path,
    )
    if not files:
        print("No files found in the bucket.")
    else:
        print(f"Files in bucket '{bucket_name}':")
        for file in files:
            print(f"- '{file}'")


@file_app.command("download")
def download_files(
    files: List[str] = typer.Argument(None, help="Space seperated list of files to download from the bucket."),
    file_list: Path | None = typer.Option(None, "--file-list", help="Text file with list of files to upload."),
    download_dir: Path = typer.Option(
        default=Path("."),
        help="Directory to download the files to.",
    ),
    bucket_version: str = typer.Option(default=None, help="Version of the bucket at which to download the files."),
    bucket_name: str = BUCKET_NAME_OPTION,
    config_path: Path = CONFIG_PATH_OPTION,
):
    """
    Download files from the bucket. This can be done by either:
        1. providing a list of files paths directly in the command line
        2. providing a directory to download the files to.
    """
    bucket_name = resolve_bucket_name(bucket_name, config_path)

    if bool(files) + bool(file_list) > 1:
        print("Please specify only one of --files or --file-list.")
        raise typer.Exit(1)

    all_files = set()
    if files:
        all_files.update(files)
    if file_list:
        with open(file_list) as f:
            for line in f:
                object_name = line.strip()
                if object_name:
                    all_files.add(object_name)

    if not all_files:
        print("No files specified for download.")
        raise typer.Exit(1)

    downloaded_files = download_files_command(
        bucket_name=bucket_name,
        all_files=list(all_files),
        download_dir=download_dir,
        bucket_version=bucket_version,
        config_path=config_path,
    )
    missing_files = all_files - set(downloaded_files)
    if missing_files:
        print("WARNING: The following files were not downloaded:")
        for file in missing_files:
            print(f"- {file}")
    else:
        print(f"The following files were downloaded to {download_dir.resolve()}:")
        for file in downloaded_files:
            print(f"- {file}")


@file_app.command("upload")
def upload_files(
    files: List[Path] | None = typer.Argument(None, help="Space seperated list of files to upload."),
    upload_dir: Path | None = typer.Option(None, "--upload-dir", help="Directory to upload all files from."),
    file_list: Path | None = typer.Option(None, "--file-list", help="Text file with list of files to upload."),
    safe_mode: Annotated[
        bool,
        typer.Option(
            "--safe-mode", help="Check if any of the files you're about to upload already exist and if so don't upload"
        ),
    ] = False,
    bucket_name: str = BUCKET_NAME_OPTION,
    config_path: Path = CONFIG_PATH_OPTION,
):
    """
    Upload files to the bucket/DivBase project by either:
        1. providing a list of files paths directly in the command line
        2. providing a directory to upload
        3. providing a text file with or a file list.
    """
    bucket_name = resolve_bucket_name(bucket_name, config_path)

    all_files = set()

    if bool(files) + bool(upload_dir) + bool(file_list) > 1:
        print("Please specify only one of --files, --upload_dir, or --file-list.")
        raise typer.Exit(1)

    if files:
        all_files.update(files)
    if upload_dir:
        all_files.update([p for p in upload_dir.iterdir() if p.is_file()])
    if file_list:
        with open(file_list) as f:
            for line in f:
                path = Path(line.strip())
                if path.is_file():
                    all_files.add(path)

    if not all_files:
        print("No files specified for upload.")
        raise typer.Exit(1)

    uploaded_files = upload_files_command(
        bucket_name=bucket_name,
        all_files=list(all_files),
        safe_mode=safe_mode,
        config_path=config_path,
    )

    if uploaded_files:
        print("Uploaded files:")
        for file in uploaded_files:
            print(f"- '{file.resolve()}'")
    else:
        print("No files were uploaded.")
