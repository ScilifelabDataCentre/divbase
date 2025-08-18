from pathlib import Path

import typer

from divbase_tools.cli_commands.config_resolver import resolve_project
from divbase_tools.cli_commands.user_config_cli import CONFIG_FILE_OPTION
from divbase_tools.s3_client import create_s3_file_manager
from divbase_tools.tasks import download_vcf_files
from divbase_tools.vcf_dimension_indexing import VCFDimensionIndexManager

PROJECT_NAME_OPTION = typer.Option(
    None,
    help="Name of the DivBase project, if not provided uses the default in your DivBase config file",
    show_default=False,
)


dimensions_app = typer.Typer(
    no_args_is_help=True,
    help="Create and inspect dimensions (number of samples, number of variants, scaffold names) of the VCF files in a project",
)

# 1. update .vcf_dimensions.yaml


@dimensions_app.command("create")
def create_dimensions_index(
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """Create a dimensions indexing file that is stored inside the project's storage bucket."""
    project_config = resolve_project(project_name=project, config_path=config_file)
    s3_file_manager = create_s3_file_manager(project_config.s3_url)
    VCFDimensionIndexManager(bucket_name=project_config.bucket_name, s3_file_manager=s3_file_manager)
    print(f"Dimensions indexing file created for project: '{project_config.name}'")

    # TODO this should check if there is a file already, and if there is it should do nothing and ask the user to run other commands instead
    # TODO submit a celery task to calculate all dimensions for VCF files


@dimensions_app.command("update")
def update_dimensions_index(
    project: str | None = PROJECT_NAME_OPTION,
    config_file: Path = CONFIG_FILE_OPTION,
):
    """Calculate and add the dimensions of a VCF file to the dimensions index file in the project."""
    project_config = resolve_project(project_name=project, config_path=config_file)
    s3_file_manager = create_s3_file_manager(project_config.s3_url)

    all_files = s3_file_manager.list_files(bucket_name=project_config.bucket_name)
    vcf_files = [file for file in all_files if file.endswith(".vcf") or file.endswith(".vcf.gz")]
    print(vcf_files)
    manager = VCFDimensionIndexManager(bucket_name=project_config.bucket_name, s3_file_manager=s3_file_manager)

    already_indexed_vcfs = manager.get_indexed_filenames()

    non_indexed_vcfs = [file for file in vcf_files if file not in already_indexed_vcfs]

    _ = download_vcf_files(
        files_to_download=non_indexed_vcfs,
        bucket_name=project_config.bucket_name,
        s3_file_manager=s3_file_manager,
    )

    for file in non_indexed_vcfs:
        manager.add_dimension_entry(vcf_filename=file)


# TODO add command to list all scaffolds available in bucket
