"""
Command line interface for managing files in a DivBase project's store on DivBase.

TODO - Download all files option.
TODO - skip checked option (aka skip files that already exist in same local dir with correct checksum).
"""

from pathlib import Path
from zoneinfo import ZoneInfo

import typer
from rich import print
from rich.table import Table
from typing_extensions import Annotated

from divbase_cli.cli_commands.version_cli import PROJECT_NAME_OPTION
from divbase_cli.cli_exceptions import UnsupportedFileNameError, UnsupportedFileTypeError
from divbase_cli.config_resolver import ensure_logged_in, resolve_download_dir, resolve_project
from divbase_cli.services.s3_files import (
    download_files_command,
    get_file_info_command,
    list_files_command,
    restore_objects_command,
    soft_delete_objects_command,
    stream_file_command,
    upload_files_command,
)
from divbase_cli.utils import format_file_size, print_rich_table_as_tsv
from divbase_lib.divbase_constants import SUPPORTED_DIVBASE_FILE_TYPES, UNSUPPORTED_CHARACTERS_IN_FILENAMES

file_app = typer.Typer(no_args_is_help=True, help="Download/upload/list files to/from the project's store on DivBase.")

NO_FILES_SPECIFIED_MSG = "No files specified for the command, exiting..."
FORMAT_AS_TSV_OPTION = typer.Option(
    False,
    "--tsv",
    help="If set, will print the output in .TSV format for easier programmatic parsing.",
)


@file_app.command("ls")
def list_files(
    format_output_as_tsv: bool = FORMAT_AS_TSV_OPTION,
    prefix_filter: str | None = typer.Option(
        None,
        "--prefix",
        "-p",
        help="Optional prefix to filter the listed files by name (only list files starting with this prefix).",
    ),
    include_results_files: bool = typer.Option(
        False,
        "--include-results-files",
        "-r",
        help="If set, will also show DivBase query results files which are hidden by default.",
    ),
    project: str | None = PROJECT_NAME_OPTION,
):
    """
    list all currently available files in the project's DivBase store.

    You can optionally filter the listed files by providing a prefix.
    By default, DivBase query results files are hidden from the listing. Use the --include-results-files option to include them.
    To see information about the versions of each file, use the 'divbase-cli files info [FILE_NAME]' command instead
    """
    project_config = resolve_project(project_name=project)
    logged_in_url = ensure_logged_in(desired_url=project_config.divbase_url)

    files = list_files_command(
        divbase_base_url=logged_in_url,
        project_name=project_config.name,
        prefix_filter=prefix_filter,
        include_results_files=include_results_files,
    )

    if not files:
        print("No files found in the project's store on DivBase.")
        return

    table = Table(title=f"Files in [bold]{project_config.name}'s [/bold] DivBase Store:")
    table.add_column("Name", justify="left", style="cyan")
    table.add_column("File size", justify="left", style="magenta", no_wrap=True)
    table.add_column("Upload date (CET)", justify="left", style="green", no_wrap=True)
    table.add_column("MD5 checksum", justify="left", style="yellow")

    for file_details in files:
        cet_timestamp = file_details.last_modified.astimezone(ZoneInfo("CET")).strftime("%Y-%m-%d %H:%M:%S %Z")
        file_size = format_file_size(size_bytes=file_details.size)
        table.add_row(
            file_details.name,
            file_size,
            cet_timestamp,
            file_details.etag,
        )

    if not format_output_as_tsv:
        print(table)
    else:
        print_rich_table_as_tsv(table=table)


@file_app.command("info")
def file_info(
    file_name: str = typer.Argument(..., help="Name of the file to get information about."),
    format_output_as_tsv: bool = FORMAT_AS_TSV_OPTION,
    project: str | None = PROJECT_NAME_OPTION,
):
    """
    Get detailed information about a specific file in the project's DivBase store.

    This includes all versions of the file and whether the file is currently marked as soft deleted.
    """
    project_config = resolve_project(project_name=project)
    logged_in_url = ensure_logged_in(desired_url=project_config.divbase_url)

    file_info = get_file_info_command(
        divbase_base_url=logged_in_url,
        project_name=project_config.name,
        object_name=file_name,
    )

    if not file_info.versions:
        # API should never return an object with no versions, but just in case
        print("No available versions for this file.")
        return

    if file_info.is_currently_deleted:
        print(
            "[bold red]Warning: This file is marked as soft deleted.\n[/bold red]"
            + "[italic]   To restore a deleted file, use the 'divbase-cli files restore' command.\n"
            + "   Or upload the file again using the 'divbase-cli files upload' command.\n[/italic]",
        )

    table = Table(
        title=f"Available Versions for '[bold]{file_info.object_name}[/bold]'",
        caption="Versions shown are ordered with the latest/current version first/at the top",
    )
    table.add_column("File size", justify="left", style="magenta", no_wrap=True)
    table.add_column("Upload date (CET)", justify="left", style="green", no_wrap=True)
    table.add_column("MD5 checksum", justify="left", style="yellow")
    table.add_column("Version ID", justify="left", style="cyan")

    for version in file_info.versions:
        cet_timestamp = version.last_modified.astimezone(ZoneInfo("CET")).strftime("%Y-%m-%d %H:%M:%S %Z")
        file_size = format_file_size(size_bytes=version.size)
        table.add_row(
            file_size,
            cet_timestamp,
            version.etag,
            version.version_id,
        )

    if not format_output_as_tsv:
        print(table)
    else:
        print_rich_table_as_tsv(table=table)


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
    project_version: str | None = typer.Option(
        default=None,
        help="User defined version of the project's at which to download the files. If not provided, downloads the latest version of all selected files.",
    ),
    project: str | None = PROJECT_NAME_OPTION,
):
    """
    Download files from the project's store on DivBase.

    This can be done by either:
        1. providing a list of files paths directly in the command line
        2. providing a text file with a list of files to download (new file on each line).

    To download the latest version of a file, just provide its name. "file1" "file2" etc.
    To download a specific/older version of a file, use the format: "file_name:version_id"
    You can get a file's version id using the 'divbase-cli file info [FILE_NAME]' command.
    You can mix and match latest and specific versions in the same command.
    E.g. to download the latest version of file1 and version "3xcdsdsdiw829x"
    of file2: 'divbase-cli files download file1 file2:3xcdsdsdiw829x'
    """
    project_config = resolve_project(project_name=project)
    logged_in_url = ensure_logged_in(desired_url=project_config.divbase_url)
    download_dir_path = resolve_download_dir(download_dir=download_dir)

    raw_files_input = _resolve_file_inputs(files=files, file_list=file_list)

    download_results = download_files_command(
        divbase_base_url=logged_in_url,
        project_name=project_config.name,
        raw_files_input=raw_files_input,
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


@file_app.command("stream")
def stream_file(
    file_name: str = typer.Argument(..., help="Name of the file you want to stream."),
    version_id: str | None = typer.Option(
        default=None,
        help="Specify this if you want to look at an older/specific version of the file. "
        "If not provided, the latest version of the file is used. "
        "To get a file's version ids, use the 'divbase-cli file info [FILE_NAME]' command.",
    ),
    project: str | None = PROJECT_NAME_OPTION,
):
    """
    Stream a file's content to standard output.

    This allows your to pipe the output to other tools like 'less', 'head', 'zcat' and 'bcftools'.

    Examples:
    - View a file: divbase-cli files stream my_file.tsv | less
    - View a gzipped file: divbase-cli files stream my_file.vcf.gz | zcat | less
    - Run a bcftools command: divbase-cli files stream my_file.vcf.gz | bcftools view -h -  # The "-" tells bcftools to read from standard input
    """
    project_config = resolve_project(project_name=project)
    logged_in_url = ensure_logged_in(desired_url=project_config.divbase_url)

    stream_file_command(
        divbase_base_url=logged_in_url,
        project_name=project_config.name,
        file_name=file_name,
        version_id=version_id,
    )


@file_app.command("upload")
def upload_files(
    files: list[Path] | None = typer.Argument(None, help="Space separated list of files to upload."),
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
):
    """
    Upload files to your project's store on DivBase:

    To provide files to upload you can either:
        1. provide a list of files paths directly in the command line
        2. provide a directory to upload
        3. provide a text file with or a file list.
    """
    project_config = resolve_project(project_name=project)
    logged_in_url = ensure_logged_in(desired_url=project_config.divbase_url)

    if bool(files) + bool(upload_dir) + bool(file_list) > 1:
        print("Please specify only one of --files, --upload_dir, or --file-list.")
        raise typer.Exit(1)

    all_files: set[Path] = set()
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
        print(NO_FILES_SPECIFIED_MSG)
        raise typer.Exit(1)

    _check_for_unsupported_files(all_files)

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


@file_app.command("rm")
def remove_files(
    files: list[str] | None = typer.Argument(
        None, help="Space seperated list of files/objects in the project's store on DivBase to delete."
    ),
    file_list: Path | None = typer.Option(None, "--file-list", help="Text file with list of files to delete."),
    dry_run: bool = typer.Option(
        False, "--dry-run", help="If set, will not actually delete the files, just print what would be deleted."
    ),
    project: str | None = PROJECT_NAME_OPTION,
):
    """
    Soft delete files from the project's store on DivBase

    To provide files to delete you can either:
        1. provide a list of file names directly in the command line
        2. provide a text file with a list of files to delete.

    Note that deleting a non existent file will be treated as a successful deletion.
    """
    project_config = resolve_project(project_name=project)
    logged_in_url = ensure_logged_in(desired_url=project_config.divbase_url)

    all_files = _resolve_file_inputs(files=files, file_list=file_list)

    if dry_run:
        print("Dry run mode enabled. The following files would have been deleted:")
        for file in all_files:
            print(f"- '{file}'")
        return

    deleted_files = soft_delete_objects_command(
        divbase_base_url=logged_in_url,
        project_name=project_config.name,
        all_files=all_files,
    )

    if deleted_files:
        print("Deleted files:")
        for file in deleted_files:
            print(f"- '{file}'")
    else:
        print("No files were deleted.")


@file_app.command("restore")
def restore_soft_deleted_files(
    files: list[str] | None = typer.Argument(
        None, help="Space seperated list of files/objects in the project's store on DivBase to restore."
    ),
    file_list: Path | None = typer.Option(None, "--file-list", help="Text file with list of files to restore."),
    project: str | None = PROJECT_NAME_OPTION,
):
    """
    Restore soft deleted files from the project's store on DivBase

    To provide files to restore you can either:
        1. provide a list of files directly in the command line.
        2. provide a text file with a list of files to restore (new file on each line).

    NOTE: Attempts to restore a file that is not soft deleted will be considered successful and the file will remain live. This means you can repeatedly run this command on the same file and get the same response.
    """
    project_config = resolve_project(project_name=project)
    logged_in_url = ensure_logged_in(desired_url=project_config.divbase_url)

    all_files = _resolve_file_inputs(files=files, file_list=file_list)

    restored_objects_response = restore_objects_command(
        divbase_base_url=logged_in_url,
        project_name=project_config.name,
        all_files=all_files,
    )

    if restored_objects_response.restored:
        print("Restored files:")
        for file in restored_objects_response.restored:
            print(f"- '{file}'")

    if restored_objects_response.not_restored:
        print("[bold red]WARNING: Some files could not be restored:[/bold red]")
        for file in restored_objects_response.not_restored:
            print(f"[red]- '{file}'[/red]")

        print(
            "Possible reasons for failed restores:\n"
            "1. The object does not exist in the bucket (e.g., a typo in the name).\n"
            "2. The object was hard-deleted and is unrecoverable.\n"
            "3. An unexpected server error occurred during the restore attempt."
        )


def _resolve_file_inputs(files: list[str] | None, file_list: Path | None) -> list[str]:
    """Helper function to resolve file inputs from command line arguments."""
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

    if not all_files:
        print(NO_FILES_SPECIFIED_MSG)
        raise typer.Exit(1)
    return list(all_files)


def _check_for_unsupported_files(all_files: set[Path]) -> None:
    """
    Helper fn to check if any of the files to be uploaded are not supported by DivBase.
    Raises error if so.

    This can be to prevent users from:
    1. Accidentally uploading unsupported files.
    2. Uploading files that have characters that we reserve for actions like filtering/querying on DivBase.
    (e.g. the syntax "[file_name]:[version_id]" can be used to download specific file versions).

    This is not a security feature, just for UX purposes.
    """
    unsupported_file_types, unsupported_chars = [], []
    for file_path in all_files:
        if not any(file_path.name.endswith(supported) for supported in SUPPORTED_DIVBASE_FILE_TYPES):
            unsupported_file_types.append(file_path)

        if any(char in file_path.name for char in UNSUPPORTED_CHARACTERS_IN_FILENAMES):
            unsupported_chars.append(file_path)

    if unsupported_file_types:
        raise UnsupportedFileTypeError(unsupported_files=unsupported_file_types)

    if unsupported_chars:
        raise UnsupportedFileNameError(unsupported_files=unsupported_chars)
