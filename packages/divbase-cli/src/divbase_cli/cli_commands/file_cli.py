"""
Command line interface for managing files in a DivBase project's storage bucket.

TODO - support for specifying versions of files when downloading files?
TODO - some duplication of logic here, but awkward as not exactly same logic for different ops.
"""

from pathlib import Path
from typing import List

import typer
from rich import print
from typing_extensions import Annotated

from divbase_cli.cli_commands.user_config_cli import CONFIG_FILE_OPTION
from divbase_cli.cli_commands.version_cli import PROJECT_NAME_OPTION
from divbase_cli.config_resolver import resolve_download_dir, resolve_project
from divbase_cli.services import (
    download_files_command,
    list_files_command,
    soft_delete_objects_command,
    upload_files_command,
)

file_app = typer.Typer(no_args_is_help=True, help="Download/upload/list files to/from the project's storage bucket.")


@file_app.command("list")
def list_files(
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """
    list all files in the project's storage bucket.

    To see files at a user specified bucket version (controlled by the 'divbase-cli version' subcommand),
    you can instead use the 'divbase version info [VERSION_NAME]' command.
    """
    project_config = resolve_project(project_name=project, config_path=config_file)
    files = list_files_command(project_config=project_config)
    if not files:
        print("No files found in the project's storage bucket.")
    else:
        print(f"Files in bucket '{project_config.bucket_name}':")
        for file in files:
            print(f"- '{file}'")


@file_app.command("download")
def download_files(
    files: List[str] = typer.Argument(None, help="Space seperated list of files/objects to download from the bucket."),
    file_list: Path | None = typer.Option(None, "--file-list", help="Text file with list of files to upload."),
    download_dir: str = typer.Option(
        None,
        help="""Directory to download the files to. 
            If not provided, defaults to what you specified in your user config. 
            If also not specified in your user config, downloads to the current directory.
            You can also specify "." to download to the current directory.""",
    ),
    bucket_version: str = typer.Option(
        default=None, help="Version of the project's storage bucket at which to download the files."
    ),
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """
    Download files from the project's storage bucket. This can be done by either:
        1. providing a list of files paths directly in the command line
        2. providing a directory to download the files to.
    """
    project_config = resolve_project(project_name=project, config_path=config_file)
    download_dir_path = resolve_download_dir(download_dir=download_dir, config_path=config_file)

    if bool(files) + bool(file_list) > 1:
        print("Please specify only one of --files or --file-list.")
        raise typer.Exit(1)

    all_files = set()
    if files:
        all_files.update(files)
    if file_list:
        with open(file_list) as f:
            for object_name in f:
                all_files.add(object_name.strip())

    if not all_files:
        print("No files specified for download.")
        raise typer.Exit(1)

    downloaded_files = download_files_command(
        project_config=project_config,
        all_files=list(all_files),
        download_dir=download_dir_path,
        bucket_version=bucket_version,
    )

    downloaded_file_names = [file.name for file in downloaded_files]
    missing_files = all_files - set(downloaded_file_names)

    if missing_files:
        print("WARNING: The following files were not downloaded:")
        for file in missing_files:
            print(f"- {file}")
    else:
        print(f"The following files were downloaded to {download_dir_path.resolve()}:")
        for file in downloaded_files:
            print(f"- '{file}'")


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
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """
    Upload files to your project's storage bucket by either:
        1. providing a list of files paths directly in the command line
        2. providing a directory to upload
        3. providing a text file with or a file list.
    """
    project_config = resolve_project(project_name=project, config_path=config_file)

    if bool(files) + bool(upload_dir) + bool(file_list) > 1:
        print("Please specify only one of --files, --upload_dir, or --file-list.")
        raise typer.Exit(1)

    all_files = set()
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
        project_config=project_config,
        all_files=list(all_files),
        safe_mode=safe_mode,
    )

    if uploaded_files:
        print("Uploaded files:")
        for object_name, file_path in uploaded_files.items():
            print(f"- '{object_name}' created from file at: '{file_path.resolve()}'")
    else:
        print("No files were uploaded.")


@file_app.command("remove")
def remove_files(
    files: List[str] | None = typer.Argument(
        None, help="Space seperated list of files/objects in the project's storage bucket to delete."
    ),
    file_list: Path | None = typer.Option(None, "--file-list", help="Text file with list of files to upload."),
    dry_run: bool = typer.Option(
        False, "--dry-run", help="If set, will not actually delete the files, just print what would be deleted."
    ),
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """
    Remove files from the project's storage bucket by either:
        1. providing a list of files paths directly in the command line
        2. providing a text file with or a file list.

    'dry_run' mode will not actually delete the files, just print what would be deleted.
    """
    project_config = resolve_project(project_name=project, config_path=config_file)

    if bool(files) + bool(file_list) > 1:
        print("Please specify only one of --files or --file-list.")
        raise typer.Exit(1)

    all_files = set()

    if files:
        all_files.update(files)
    if file_list:
        with open(file_list) as f:
            for line in f:
                all_files.add(line.strip())

    if dry_run:
        print("Dry run mode enabled. The following files would have been deleted:")
        for file in all_files:
            print(f"- '{file}'")
        return

    deleted_files = soft_delete_objects_command(
        project_config=project_config,
        all_files=list(all_files),
    )

    if deleted_files:
        print("Deleted files:")
        for file in deleted_files:
            print(f"- '{file}'")
    else:
        print("No files were deleted.")
