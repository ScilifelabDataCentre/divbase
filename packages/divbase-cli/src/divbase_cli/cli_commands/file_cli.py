"""
Command line interface for managing files in a DivBase project's store on DivBase.

TODO - support for specifying versions of files when downloading files?
TODO - Download all files option.
TODO - skip checked option (aka skip files that already exist in same local dir with correct checksum).
"""

from pathlib import Path

import typer
from rich import print
from typing_extensions import Annotated

from divbase_cli.cli_commands.user_config_cli import CONFIG_FILE_OPTION
from divbase_cli.cli_commands.version_cli import PROJECT_NAME_OPTION
from divbase_cli.config_resolver import ensure_logged_in, resolve_download_dir, resolve_project
from divbase_cli.services import (
    download_files_command,
    list_files_command,
    soft_delete_objects_command,
    upload_files_command,
)

file_app = typer.Typer(no_args_is_help=True, help="Download/upload/list files to/from the project's store on DivBase.")


@file_app.command("list")
def list_files(
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """
    list all files in the project's DivBase store.

    To see files at a user specified project version (controlled by the 'divbase-cli version' subcommand),
    you can instead use the 'divbase-cli version info [VERSION_NAME]' command.
    """
    project_config = resolve_project(project_name=project, config_path=config_file)
    logged_in_url = ensure_logged_in(config_path=config_file, desired_url=project_config.divbase_url)

    files = list_files_command(divbase_base_url=logged_in_url, project_name=project_config.name)
    if not files:
        print("No files found in the project's store on DivBase.")
    else:
        print(f"Files in the project '{project_config.name}':")
        for file in files:
            print(f"- '{file}'")


@file_app.command("download")
def download_files(
    files: list[str] = typer.Argument(
        None, help="Space separated list of files/objects to download from the project's store on DivBase."
    ),
    file_list: Path | None = typer.Option(None, "--file-list", help="Text file with list of files to upload."),
    download_dir: str = typer.Option(
        None,
        help="""Directory to download the files to. 
            If not provided, defaults to what you specified in your user config. 
            If also not specified in your user config, downloads to the current directory.
            You can also specify "." to download to the current directory.""",
    ),
    disable_verify_checksums: Annotated[
        bool,
        typer.Option(
            "--disable-verify-checksums",
            help="Turn off checksum verification which is on by default. "
            "Checksum verification means all downloaded files are verified against their MD5 checksums."
            "It is recommended to leave checksum verification enabled unless you have a specific reason to disable it.",
        ),
    ] = False,
    project_version: str = typer.Option(
        default=None,
        help="User defined version of the project's at which to download the files. If not provided, downloads the latest version of all selected files.",
    ),
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """
    Download files from the project's store on DivBase. This can be done by either:
        1. providing a list of files paths directly in the command line
        2. providing a directory to download the files to.
    """
    project_config = resolve_project(project_name=project, config_path=config_file)
    logged_in_url = ensure_logged_in(config_path=config_file, desired_url=project_config.divbase_url)
    download_dir_path = resolve_download_dir(download_dir=download_dir, config_path=config_file)

    if bool(files) + bool(file_list) > 1:
        print("Please specify only one of --files or --file-list.")
        raise typer.Exit(1)

    all_files: set[str] = set()
    if files:
        all_files.update(files)
    if file_list:
        with open(file_list) as f:
            for object_name in f:
                all_files.add(object_name.strip())

    if not all_files:
        print("No files specified for download.")
        raise typer.Exit(1)

    download_results = download_files_command(
        divbase_base_url=logged_in_url,
        project_name=project_config.name,
        all_files=list(all_files),
        download_dir=download_dir_path,
        verify_checksums=not disable_verify_checksums,
        project_version=project_version,
    )

    if download_results.successful:
        print("[green bold]Successfully downloaded the following files:[/green bold]")
        for success in download_results.successful:
            print(f"- '{success.object_name}' downloaded to: '{success.file_path.resolve()}'")
    if download_results.failed:
        print("[red bold]ERROR: Failed to download the following files:[/red bold]")
        for failed in download_results.failed:
            print(f"[red]- '{failed.object_name}': Exception: '{failed.exception}'[/red]")

        raise typer.Exit(1)


@file_app.command("upload")
def upload_files(
    files: list[Path] | None = typer.Argument(None, help="Space seperated list of files to upload."),
    upload_dir: Path | None = typer.Option(None, "--upload-dir", help="Directory to upload all files from."),
    file_list: Path | None = typer.Option(None, "--file-list", help="Text file with list of files to upload."),
    disable_safe_mode: Annotated[
        bool,
        typer.Option(
            "--disable-safe-mode",
            help="Turn off safe mode which is on by default. Safe mode adds 2 extra bits of security by first calculating the MD5 checksum of each file that you're about to upload:"
            "(1) Checks if any of the files you're about to upload already exist (by comparing name and checksum) and if so stops the upload process."
            "(2) Sends the file's checksum when the file is uploaded so the server can verify the upload was successful (by calculating and comparing the checksums)."
            "It is recommended to leave safe mode enabled unless you have a specific reason to disable it.",
        ),
    ] = False,
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """
    Upload files to your project's store on DivBase by either:
        1. providing a list of files paths directly in the command line
        2. providing a directory to upload
        3. providing a text file with or a file list.
    """
    project_config = resolve_project(project_name=project, config_path=config_file)
    logged_in_url = ensure_logged_in(config_path=config_file, desired_url=project_config.divbase_url)

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

    uploaded_results = upload_files_command(
        project_name=project_config.name,
        divbase_base_url=logged_in_url,
        all_files=list(all_files),
        safe_mode=not disable_safe_mode,
    )

    if uploaded_results.successful:
        print("[green bold] The following files were successfully uploaded: [/green bold]")
        for object in uploaded_results.successful:
            print(f"- '{object.object_name}' created from file at: '{object.file_path.resolve()}'")

    if uploaded_results.failed:
        print("[red bold]ERROR: Failed to upload the following files:[/red bold]")
        for failed in uploaded_results.failed:
            print(f"[red]- '{failed.object_name}': Exception: '{failed.exception}'[/red]")

        raise typer.Exit(1)


@file_app.command("remove")
def remove_files(
    files: list[str] | None = typer.Argument(
        None, help="Space seperated list of files/objects in the project's store on DivBase to delete."
    ),
    file_list: Path | None = typer.Option(None, "--file-list", help="Text file with list of files to upload."),
    dry_run: bool = typer.Option(
        False, "--dry-run", help="If set, will not actually delete the files, just print what would be deleted."
    ),
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """
    Remove files from the project's store on DivBase by either:
        1. providing a list of files paths directly in the command line
        2. providing a text file with or a file list.

    'dry_run' mode will not actually delete the files, just print what would be deleted.
    """
    project_config = resolve_project(project_name=project, config_path=config_file)
    logged_in_url = ensure_logged_in(config_path=config_file, desired_url=project_config.divbase_url)

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
        divbase_base_url=logged_in_url,
        project_name=project_config.name,
        all_files=list(all_files),
    )

    if deleted_files:
        print("Deleted files:")
        for file in deleted_files:
            print(f"- '{file}'")
    else:
        print("No files were deleted.")
