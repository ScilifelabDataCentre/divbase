"""
Command line interface for managing files in a DivBase project's store on DivBase.
"""

from glob import glob
from pathlib import Path
from zoneinfo import ZoneInfo

import typer
from rich import print
from rich.console import Console
from rich.table import Table
from typing_extensions import Annotated

from divbase_cli.cli_commands.shared_args_options import DOWNLOAD_DIR_OPTION, FORMAT_AS_TSV_OPTION, PROJECT_NAME_OPTION
from divbase_cli.config_resolver import ensure_logged_in, resolve_download_dir, resolve_project
from divbase_cli.services.project_versions import get_version_details_command
from divbase_cli.services.s3_files import (
    ToDownload,
    ToUpload,
    download_files_command,
    filter_out_already_downloaded_files,
    get_file_info_command,
    list_files_command,
    list_soft_deleted_files_command,
    make_directories_command,
    restore_objects_command,
    soft_delete_objects_command,
    stream_file_command,
    upload_files_command,
)
from divbase_cli.utils import print_rich_table_as_tsv
from divbase_lib.divbase_constants import (
    SUPPORTED_DIVBASE_FILE_TYPES,
    UNSUPPORTED_CHARACTERS_DISPLAY,
    UNSUPPORTED_CHARACTERS_IN_FILENAMES,
)
from divbase_lib.utils import format_file_size

file_app = typer.Typer(no_args_is_help=True, help="Download/upload/list files to/from the project's store on DivBase.")

NO_FILES_SPECIFIED_MSG = "No files specified for the command, exiting..."
NO_UPLOAD_MATCHES_MSG = "Error: The following file paths or glob patterns did not match any existing files:"

DISABLE_VERIFY_CHECKSUMS_OPTION = typer.Option(
    False,
    "--disable-verify-checksums",
    "-nc",
    help="Turn off checksum verification which is on by default. "
    "Checksum verification means all downloaded files are verified against their MD5 checksums. "
    "It is recommended to leave checksum verification enabled unless you have a specific reason to disable it.",
)
PROJECT_VERSION_OPTION = typer.Option(
    None,
    "--project-version",
    "-pv",
    help="User defined version of the project's at which to download the files. If not provided, downloads the latest version of all selected files.",
)
DRY_RUN_OPTION = typer.Option(
    False,
    "--dry-run",
    "-n",
    help="If set, will not actually download the files, just print what would be downloaded.",
)


@file_app.command("ls")
def list_files(
    prefix: str | None = typer.Argument(
        None,
        help="Filter by prefix, if you want to see files and folders inside a specific dir, include the '/' at the end: e.g. 'vcfs/' . ",
    ),
    detailed: bool = typer.Option(
        False,
        "--detailed",
        "-l",
        help="Show a detailed view including the file size and upload date.",
    ),
    format_output_as_tsv: bool = FORMAT_AS_TSV_OPTION,
    include_results_files: bool = typer.Option(
        False,
        "--include-results-files",
        "-r",
        help="If set, will also show DivBase query results files which are hidden by default.",
    ),
    show_deleted_files: bool = typer.Option(
        False,
        "--show-deleted-files",
        "-s",
        help="Show the files in the project that are currently soft deleted. These files can be recovered within a certain time frame after deletion",
    ),
    project: str | None = PROJECT_NAME_OPTION,
):
    """
    List files and folders in the project's DivBase store.
    # TODO - support a recursive option to show all files and subfolders as well? - can be called --tree or --recursive?

    Examples:
    - List all files and folders in the project:
        divbase-cli files ls
    - List files and folders in the 'vcfs/' folder:
        divbase-cli files ls vcfs/
    - List all files and folders starting with 'sample':
        divbase-cli files ls sample
    - List all files including DivBase query results files (hidden by default):
        divbase-cli files ls --include-results-files
    """
    if show_deleted_files and include_results_files:
        raise typer.BadParameter(
            message="The --show-deleted-files option cannot be used with the flag --include-results-files. Please use these options separately."
        )

    project_config = resolve_project(project_name=project)
    logged_in_url = ensure_logged_in(desired_url=project_config.divbase_url)

    if show_deleted_files:
        files = list_soft_deleted_files_command(
            divbase_base_url=logged_in_url,
            project_name=project_config.name,
            prefix=prefix,
        )

        if not files:
            print(f"No soft deleted files found for the project '{project_config.name}'.")
            return

        print(f"Soft deleted files for the project '{project_config.name}':")
        for file_details in files:
            cet_timestamp = file_details.last_modified.astimezone(ZoneInfo("CET")).strftime("%Y-%m-%d %H:%M:%S %Z")
            if file_details.name.endswith("/"):
                print(f"- [bold blue]'{file_details.name}'[/bold blue] (deleted at: '{cet_timestamp}')")
            else:
                print(f"- '{file_details.name}' (deleted at: '{cet_timestamp}')")
        print("\nTo restore a soft deleted file, use the 'divbase-cli files restore' command.")
        print("To get more information about one of the soft deleted files, use the 'divbase-cli files info' command.")
        return

    files, folders = list_files_command(
        divbase_base_url=logged_in_url,
        project_name=project_config.name,
        prefix=prefix,
        include_results_files=include_results_files,
    )

    if not files and not folders:
        print("No files or folders found in the project's store on DivBase.")
        return

    if detailed:
        _print_ls_detailed(files=files, folders=folders, format_as_tsv=format_output_as_tsv)
    else:
        _print_ls_simple(files=files, folders=folders)


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
    file_list: Path | None = typer.Option(None, "--file-list", help="Text file with list of files to download."),
    download_dir: str = DOWNLOAD_DIR_OPTION,
    dry_run: bool = DRY_RUN_OPTION,
    disable_verify_checksums: bool = DISABLE_VERIFY_CHECKSUMS_OPTION,
    project_version: str | None = PROJECT_VERSION_OPTION,
    project: str | None = PROJECT_NAME_OPTION,
):
    """
    Download files from the project's store on DivBase.

    This can be done by either:
        1. providing a list of files paths directly in the command line
        2. providing a text file with a list of files to download (new file on each line).

    To download the latest version of a file, just provide its name. "file1" "file2" etc.
    To download a specific/older version of a file, use the format: "file_name:version_id"
    You can get a file's version id using the 'divbase-cli files info [FILE_NAME]' command.
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
        dry_run=dry_run,
        project_version=project_version,
    )

    _pretty_print_download_results(download_results=download_results)


@file_app.command("download-all")
def download_all_files(
    download_dir: str = DOWNLOAD_DIR_OPTION,
    resume: bool = typer.Option(
        False,
        "--resume",
        "-r",
        help="If set, will attempt to resume an interrupted download. Will check which files have already been fully downloaded (by checking if a file with the same name and checksum already exists in the download directory) and skip downloading those files again.",
    ),
    dry_run: bool = DRY_RUN_OPTION,
    disable_verify_checksums: bool = DISABLE_VERIFY_CHECKSUMS_OPTION,
    project_version: str | None = PROJECT_VERSION_OPTION,
    project: str | None = PROJECT_NAME_OPTION,
):
    """
    Download all files in the project's store on DivBase.
    Before the download proceeds you'll be prompted if you want to continue.
    DivBase Query results files will not be included in the download.

    You can resume ('--resume' / '-r') a 'download-all' command, just make sure you're downloading into the same directory.
    """
    if resume and disable_verify_checksums:
        print(
            "The --resume and --disable-verify-checksums options cannot be used together, "
            "as checksums are used to determine which files don't need to be downloaded. \n"
            "Exiting... "
        )
        return

    project_config = resolve_project(project_name=project)
    logged_in_url = ensure_logged_in(desired_url=project_config.divbase_url)
    download_dir_path = resolve_download_dir(download_dir=download_dir)

    all_files: list[ToDownload] = []
    if project_version:
        version_details = get_version_details_command(
            project_name=project_config.name, divbase_base_url=logged_in_url, version_name=project_version
        )
        for version_name, file_details in version_details.files.items():
            all_files.append(
                ToDownload(
                    name=version_name,
                    etag=file_details["etag"],
                    size_bytes=file_details["size"],
                    version_id=file_details["version_id"],
                )
            )

    else:
        files, _ = list_files_command(
            divbase_base_url=logged_in_url,
            project_name=project_config.name,
            prefix=None,
            include_results_files=False,
            file_system_view=False,  # we want all files in a flat list to download al of them, not a file system view with folders and files
        )
        for file_details in files:
            all_files.append(
                ToDownload(
                    name=file_details.name,
                    etag=file_details.etag,
                    size_bytes=file_details.size,
                    version_id=None,  # latest version
                )
            )

    if not all_files:
        print("No files to download as there are no files in the project's store.")
        return

    # filter files to download based on those which already exist.
    if resume:
        files_to_download, files_to_overwrite = filter_out_already_downloaded_files(
            all_files=all_files, download_dir=download_dir_path
        )
        if files_to_overwrite:
            print(
                "[yellow bold]Warning: The following files already exist in the download directory but have a different checksum."
                "If you choose to proceed, these files will be overwritten by the download:[/yellow bold]"
            )
            for file in files_to_overwrite:
                print(f"- '{file.name}'")

        files_to_download.extend(files_to_overwrite)

        if not files_to_download:
            print("No files left to download, your folder matches the project's store on DivBase, exiting...")
            return
    else:
        files_to_download = all_files

    total_size_bytes = sum(file.size_bytes for file in files_to_download)
    formatted_total_size = format_file_size(size_bytes=total_size_bytes)

    print(f"There are '{len(files_to_download)}' files to download with a total size of: {formatted_total_size}.")

    if not dry_run:
        # (dry run will auto exit in the download_files_command)
        do_download = typer.confirm("Do you want to proceed with the download?")
        if not do_download:
            print("Download cancelled...")
            return

    if project_version:
        raw_files_input = [f"{file.name}:{file.version_id}" for file in files_to_download]
    else:
        raw_files_input = [file.name for file in files_to_download]
    download_results = download_files_command(
        divbase_base_url=logged_in_url,
        project_name=project_config.name,
        raw_files_input=raw_files_input,
        download_dir=download_dir_path,
        verify_checksums=not disable_verify_checksums,
        dry_run=dry_run,
        project_version=None,  # We already know the version id of each file, so can skip the processing here.
    )

    _pretty_print_download_results(download_results=download_results)


@file_app.command("stream")
def stream_file(
    file_name: str = typer.Argument(..., help="Name of the file you want to stream."),
    version_id: str | None = typer.Option(
        default=None,
        help="Specify this if you want to look at an older/specific version of the file. "
        "If not provided, the latest version of the file is used. "
        "To get a file's version ids, use the 'divbase-cli files info [FILE_NAME]' command.",
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
    files: list[str] | None = typer.Argument(None, help="Space separated list of files or glob patterns to upload."),
    file_list: Path | None = typer.Option(None, "--file-list", "-l", help="Text file with list of files to upload."),
    remote_dir: str | None = typer.Option(
        None,
        "--to",
        "-t",
        help="Remote folder to upload into, e.g. 'vcfs/batch1/'. Folder does not need to exist before uploading.",
    ),
    skip_existing: bool = typer.Option(
        False,
        "--skip-existing",
        "-s",
        help="If set, will skip already uploaded files, from the files you provided in the command. "
        "Already uploaded files are determined by checking if a file with the same name and MD5 checksum already exists in the project's store on DivBase. ",
    ),
    recursive: bool = typer.Option(
        False,
        "--recursive",
        "-r",
        "-R",
        help="If set, recursively include subdirectories contents when uploading (i.e. '**' is expanded). "
        "Without this flag, patterns only match files in the specified directory. "
        "Put your argument in quotes to prevent shell expansion of the glob before it gets to the CLI, e.g. files upload 'data/**/*.vcf.gz' ",
    ),
    dry_run: bool = typer.Option(
        False, "--dry-run", "-n", help="If set, will show what files would be uploaded, but not actually upload them."
    ),
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
    Upload files to your project's store on DivBase.

    By default only the file name is used when storing the file in the project's store. Use '--to' to place files inside a remote folder.
    To upload files recursively, see the --recursive flag.

    Examples:
        # Upload multiple files by specifying them one after another
        divbase-cli files upload file1.vcf.gz file2.tsv

        # Upload all .vcf.gz files in current directory into a remote directory called 'experiment1/'
        divbase-cli files upload "*.vcf.gz" --to experiment1/

        # Upload all files in a directory (non-recursive) to the root of the project store
        divbase-cli files upload "/path/to/data/*"

        # Upload all files inside a directory and its subdirectories (recursive) into a remote directory called 'experiment1/'
        divbase-cli files upload --recursive "/path/to/data/**" --to experiment1/

        # Upload from a text file list (one file path per line) into a remote directory called 'experiment1/'
        divbase-cli files upload --file-list files_to_upload.txt --to experiment1/
    """
    project_config = resolve_project(project_name=project)
    logged_in_url = ensure_logged_in(desired_url=project_config.divbase_url)

    if bool(files) + bool(file_list) > 1:
        print("You cannot specify files as arguments and provide a --file-list.")
        raise typer.Exit(1)

    if skip_existing and disable_safe_mode:
        print(
            "The --skip-existing and --disable-safe-mode options cannot be used together.\n"
            "Safe mode calculates file checksums, "
            "and these checksums are needed for --skip-existing to know what files to skip uploading."
        )
        raise typer.Exit(1)

    to_upload = _resolve_user_upload_inputs(
        files=files,
        file_list=file_list,
        recursive=recursive,
        remote_dir=remote_dir,
    )
    _check_for_duplicate_object_keys(to_upload)
    _check_for_unsupported_files(to_upload)

    uploaded_results = upload_files_command(
        project_name=project_config.name,
        divbase_base_url=logged_in_url,
        all_files=to_upload,
        safe_mode=not disable_safe_mode,
        skip_existing=skip_existing,
        dry_run=dry_run,
    )

    if uploaded_results.skipped:
        print("[yellow bold]\nThe following files were skipped: [/yellow bold]")
        for file in uploaded_results.skipped:
            print(f"[yellow]- '{file.object_name}' reason: {file.reason}[/yellow]")

    if uploaded_results.successful:
        print("[green bold]\nThe following files were successfully uploaded: [/green bold]")
        for file in uploaded_results.successful:
            print(f"[green]- '{file.object_name}' created from file at: '{file.file_path.resolve()}'[/green]")

    if uploaded_results.failed:
        print("[red bold]\nERROR: Failed to upload the following files:[/red bold]")
        for failed in uploaded_results.failed:
            print(f"[red]- '{failed.object_name}': Exception: '{failed.exception}'[/red]")
        raise typer.Exit(1)
    elif not uploaded_results.successful:
        print("\nNo files were uploaded.")
    else:
        print("\n[green bold]All files uploaded successfully![/green bold]")


@file_app.command("mkdir")
def make_directory(
    directories: list[str] = typer.Argument(
        ..., help="space separated list of directories to create in the project's store on DivBase."
    ),
    project: str | None = PROJECT_NAME_OPTION,
):
    """
    Create directory(ies) in your project store.

    Provide a single directory name or multiple directory names separated by space. E.g. 'dir1' or 'dir1 dir2 dir3'.
    """
    project_config = resolve_project(project_name=project)
    logged_in_url = ensure_logged_in(desired_url=project_config.divbase_url)

    cleaned_dir_names = _sanitize_directory_names(directories)
    dirs_created = make_directories_command(
        directories=cleaned_dir_names,
        divbase_base_url=logged_in_url,
        project_name=project_config.name,
    )
    if dirs_created.failed:
        print("[red bold]WARNING: Failed to create the following directories:[/red bold]")
        for dir in dirs_created.failed:
            print(f"[red]'{dir}'[/red]")
        print("This could be due to invalid characters in the directory name or an unexpected server error.\n")

    if dirs_created.created:
        print("Successfully created the following directories:")
        for dir in dirs_created.created:
            print(f"[bold blue]'{dir}'[/bold blue]")


@file_app.command("rmdir")
def remove_directory(
    directories: list[str] = typer.Argument(
        ..., help="space separated list of directories to remove from the project's store on DivBase."
    ),
    project: str | None = PROJECT_NAME_OPTION,
):
    """
    Remove an empty directory from your project store.

    NOTE:
    - Any files inside the directory must be removed before you can remove the directory.
    - If the directory does not exist, the command will be treated as a successful deletion.
    """
    project_config = resolve_project(project_name=project)
    logged_in_url = ensure_logged_in(desired_url=project_config.divbase_url)

    directories = _sanitize_directory_names(directories)
    for dir in directories:
        files, _ = list_files_command(
            divbase_base_url=logged_in_url,
            project_name=project_config.name,
            prefix=dir,
            include_results_files=True,
            file_system_view=False,
        )
        non_placeholder_files = [f for f in files if f.name != dir]
        if non_placeholder_files:
            print(
                f"[red bold]ERROR: The directory '{dir}' is not empty. \n"
                f"You must first remove all files inside it first using 'divbase-cli files rm'.[/red bold]"
            )
            raise typer.Exit(1)

    # NOTE: as dirs are just empty files which end with "/", we can just do a regular soft delete on them.
    deleted_dirs = soft_delete_objects_command(
        divbase_base_url=logged_in_url,
        project_name=project_config.name,
        all_files=directories,
    )

    if deleted_dirs:
        print("Deleted the following directories:")
        for dir in deleted_dirs:
            print(f"[bold blue]'{dir}'[/bold blue]")


@file_app.command("rm")
def remove_files(
    files: list[str] | None = typer.Argument(
        None, help="Space separated list of files/objects in the project's store on DivBase to delete."
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
        None, help="Space separated list of files/objects in the project's store on DivBase to restore."
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


def _print_ls_simple(files, folders) -> None:
    """Compact folder-aware listing: bold blue folders first, then plain filenames."""
    console = Console()
    for folder in folders:
        console.print(f"[bold blue]{folder}[/bold blue]")
    for file_details in files:
        console.print(file_details.name)


def _print_ls_detailed(files, folders, format_as_tsv: bool) -> None:
    """Unix ls -l style: one entry per line with size, date, checksum, and name."""
    console = Console()

    if format_as_tsv:
        for folder in folders:
            folder_name = folder[:-1] if folder.endswith("/") else folder
            console.print(f"{folder_name}\t-\t-\t")
        for f in files:
            cet = f.last_modified.astimezone(ZoneInfo("CET")).strftime("%Y-%m-%d %H:%M:%S %Z")
            console.print(f"{f.name}\t{format_file_size(f.size)}\t{cet}")
        return

    for folder in folders:
        console.print(f"[bold blue]{folder}[/bold blue]")
    for f in files:
        size_str = format_file_size(f.size).rjust(10)
        cet = f.last_modified.astimezone(ZoneInfo("CET")).strftime("%Y-%m-%d %H:%M:%S %Z")
        console.print(f"{f.name} [magenta]{size_str}[/magenta]  [green]{cet}[/green]")


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


def _resolve_user_upload_inputs(
    files: list[str] | None,
    file_list: Path | None,
    recursive: bool,
    remote_dir: str | None,
) -> list[ToUpload]:
    """
    Take user CLI input arguments for files upload cmd and create a list of ToUpload objects.

    Handles the multiple ways a user can specify files for upload.
    """
    to_upload: list[ToUpload] = []
    missing_patterns: list[str] = []
    seen: set[Path] = set()

    if remote_dir:
        remote_dir = remote_dir[1:] if remote_dir.startswith("/") else remote_dir
        remote_dir = remote_dir if remote_dir.endswith("/") else remote_dir + "/"
        if any(char in remote_dir for char in UNSUPPORTED_CHARACTERS_IN_FILENAMES):
            print(
                f"[red bold]Error: The --to destination path contains unsupported characters.[/red bold]\n"
                f"[red]Unsupported characters: {UNSUPPORTED_CHARACTERS_DISPLAY}[/red]\n"
            )
            raise typer.Exit(1)

    if recursive and files and not any("*" in p for p in files) and any(len(Path(p).parts) > 1 for p in files):
        print(
            "[red bold]Error: --recursive was passed but all arguments provided are file paths.[/red bold]\n"
            "Your shell likely expanded the glob before passing it to this command.\n"
            "To solve the problem, wrap the pattern(s) in quotes to prevent shell expansion, e.g.:\n"
            "  [green]divbase-cli files upload 'my_data/**/*.vcf.gz' --recursive[/green] AND NOT \n"
            "  [red]divbase-cli files upload my_data/**/*.vcf.gz --recursive[/red]"
        )
        raise typer.Exit(1)

    if files:
        for pattern in files:
            matched = [Path(p) for p in glob(pattern, recursive=recursive) if Path(p).is_file()]
            if not matched:
                if Path(pattern).is_dir():
                    continue  # shell-expanded a directory into the args list; skip it
                missing_patterns.append(pattern)
                continue
            for file in matched:
                if file in seen:
                    continue
                seen.add(file)
                if recursive and "**" in pattern:
                    root_str = pattern.split("**")[0]
                    root = Path(root_str) if root_str else Path(".")
                    key = str(file.relative_to(root)).replace("\\", "/")
                else:
                    key = file.name

                destination_key = (remote_dir or "") + key
                to_upload.append(ToUpload(file_path=file, destination_key=destination_key))

    if file_list:
        missing_files: list[Path] = []
        with open(file_list) as f:
            for line in f:
                if not line.strip():
                    continue
                path = Path(line.strip())
                if not path.is_file():
                    missing_files.append(path)
                    continue
                if path not in seen:
                    seen.add(path)
                    to_upload.append(ToUpload(file_path=path, destination_key=(remote_dir or "") + path.name))

        if missing_files:
            print("[red bold]Error: The following file paths provided in your --file-list were not found:[/red bold]")
            for path in missing_files:
                print(f"[red]- '{path}'[/red]")
            raise typer.Exit(1)

    if missing_patterns:
        print(f"[red bold]{NO_UPLOAD_MATCHES_MSG}[/red bold]")
        for pattern in missing_patterns:
            print(f"[red]- '{pattern}'[/red]")
        raise typer.Exit(1)

    if not to_upload:
        print(NO_FILES_SPECIFIED_MSG)
        raise typer.Exit(1)

    return to_upload


def _check_for_duplicate_object_keys(to_upload: list[ToUpload]) -> None:
    """Check that no two files resolve to the same S3 object key — they would silently overwrite each other."""
    key_to_paths: dict[str, list[Path]] = {}
    for s in to_upload:
        key_to_paths.setdefault(s.destination_key, []).append(s.file_path)

    duplicates = {key: paths for key, paths in key_to_paths.items() if len(paths) > 1}
    if duplicates:
        print(
            "[red bold]Error: The following S3 destination keys appear more than once in the upload list.[/red bold]\n"
            "[red]Multiple local files would overwrite each other at the same key:[/red]\n"
        )
        for key, paths in duplicates.items():
            print(f"[red bold]'{key}':[/red bold]")
            for p in paths:
                print(f"[red]- {p.resolve()}[/red]")
        print("\nEnsure each file maps to a unique destination key, or use --to to place them in separate folders.")
        raise typer.Exit(1)


def _check_for_unsupported_files(all_files: list[ToUpload]) -> None:
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
    for file in all_files:
        if not any(file.file_name.endswith(supported) for supported in SUPPORTED_DIVBASE_FILE_TYPES):
            unsupported_file_types.append(file.file_name)

        if any(char in file.destination_key for char in UNSUPPORTED_CHARACTERS_IN_FILENAMES):
            unsupported_chars.append(file.file_name)

    if unsupported_file_types:
        print(
            f"[red bold]Error: The following files' types are not supported by DivBase and therefore cannot be uploaded: [/red bold] \n"
            f"[red]{'\n'.join(str(file) for file in unsupported_file_types)}\n\n[/red]"
            f"DivBase currently supports the following file types: '{', '.join(SUPPORTED_DIVBASE_FILE_TYPES)}'\n"
            "Note that while you can upload '.tbi' and '.csi' files they are not used by DivBase in queries, we create our own index files instead. \n"
            "If you want us to support another file type, please let us know.",
        )

    if unsupported_chars:
        print(
            f"[red bold]Error: The following file(s) cannot be uploaded because their name or destination path contains unsupported characters:[/red bold]\n"
            f"[red]{'\n'.join(str(file) for file in unsupported_chars)}\n\n[/red]"
            f"Unsupported characters: [red]{UNSUPPORTED_CHARACTERS_DISPLAY}[/red]\n"
            "Please rename the file or choose a different --to destination and try again.",
        )

    if unsupported_file_types or unsupported_chars:
        raise typer.Exit(1)


def _sanitize_directory_names(directories: list[str]) -> list[str]:
    cleaned_dir_names = []
    for directory in directories:
        if any(char in directory for char in UNSUPPORTED_CHARACTERS_IN_FILENAMES):
            print(
                f"[red bold]Error: The directory name '{directory}' contains unsupported characters.[/red bold]\n"
                f"[red]Unsupported characters: {UNSUPPORTED_CHARACTERS_DISPLAY}[/red]\n"
            )
            raise typer.Exit(1)

        directory = directory if directory.endswith("/") else directory + "/"
        directory = directory[1:] if directory.startswith("/") else directory
        cleaned_dir_names.append(directory)
    return cleaned_dir_names


def _pretty_print_download_results(download_results):
    """Helper fn used by download and download all commands to print the results of the download in a nice format."""
    if download_results.successful:
        print("\n[green bold]Successfully downloaded the following files:[/green bold]")
        for success in download_results.successful:
            print(f"- '{success.object_name}' downloaded to: '{success.file_path.resolve()}'")

    if download_results.failed:
        print("\n[red bold]ERROR: Failed to download the following files:[/red bold]")
        for failed in download_results.failed:
            print(f"[red]- '{failed.object_name}': [/red] Exception: '{failed.exception}'")
        raise typer.Exit(1)
