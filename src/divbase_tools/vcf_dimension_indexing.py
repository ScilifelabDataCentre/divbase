import datetime
import gzip
import logging
import os
from dataclasses import dataclass, field
from pathlib import Path

import botocore
import typer
import yaml

from divbase_tools.exceptions import ObjectDoesNotExistError
from divbase_tools.s3_client import S3FileManager, create_s3_file_manager
from divbase_tools.user_config import ProjectConfig

DIMENSIONS_FILE_NAME = ".vcf_dimensions.yaml"
DEFAULT_CONFIG_PATH = Path.home() / ".config" / ".divbase_tools.yaml"

CONFIG_FILE_OPTION = typer.Option(
    DEFAULT_CONFIG_PATH,
    "--config",
    "-c",
    help="Path to your user configuration file. By default it is stored at ~/.config/.divbase_tools.yaml.",
)

PROJECT_NAME_OPTION = typer.Option(
    None,
    help="Name of the DivBase project, if not provided uses the default in your DivBase config file",
    show_default=False,
)
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


os.environ["DIVBASE_ENV"] = "local"  ## to ensure connection to local for now


@dataclass
class VCFDimensionIndexManager:
    bucket_name: str
    s3_file_manager: S3FileManager

    # not defined till post_init method run.
    dimensions_info: dict = field(init=False)

    def __post_init__(self):
        """
        Creates an empty .vcf_dimensions.yaml file in the bucket upon object init.
        """
        self.dimensions_info = self._get_bucket_dimensions_file()

        if self.dimensions_info:
            logger.info(f"VCF_dimensions index file found inside Bucket: {self.bucket_name}.")
        else:
            logger.info(f"Creating a new VCF_dimensions index file in bucket: {self.bucket_name}.")
            yaml_data = {"dimensions": []}
            self._upload_bucket_dimensions_file(version_data=yaml_data)
            self.dimensions_info = yaml_data

    def add_dimension_entry(self, vcf_filename: str) -> None:
        """
        Append a new dimension entry to .vcf_dimensions.yaml if not already present for that VCF file.
        Calls submethods to calculate dimensions and upload the updated YAML file to bucket.
        """
        yaml_data = self.dimensions_info
        dimensions = yaml_data.get("dimensions", [])

        if any(entry.get("filename") == vcf_filename for entry in dimensions):
            print(f"Filename '{vcf_filename}' is already present in .vcf_dimensions.yaml.")
            return

        vcf_path = Path(vcf_filename)
        vcf_dims = self._wrapper_calculate_dimensions(vcf_path)

        timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
        new_entry = {
            "filename": vcf_filename,
            "timestamp": timestamp,
            # "version": "TODO: THE BUCKET VERSION CHECKSUM SHOULD BE HERE. THEN, ONE RECORD PER VERSION OF EACH FILE",
            "dimensions": vcf_dims,
        }
        dimensions.append(new_entry)
        yaml_data["dimensions"] = dimensions

        self._upload_bucket_dimensions_file(version_data=yaml_data)
        print(f"Added new entry for {vcf_filename} to .vcf_dimensions.yaml.")

    def get_indexed_filenames(self) -> list[str]:
        """
        Returns a list of all filenames already indexed in .vcf_dimensions.yaml for this bucket.
        """
        yaml_data = self._get_bucket_dimensions_file()
        return [entry.get("filename") for entry in yaml_data.get("dimensions", []) if "filename" in entry]

    def _wrapper_calculate_dimensions(self, vcf_path: Path) -> dict:
        """
        Reads compressed or uncompressed VCF files and returns their dimensions.

        Two different entry points to self._extract_dimensions_from_opened_vcf() in
        order to comply with Ruff linting of context managers.
        """
        print(f"Reading: {vcf_path} ...")
        try:
            with gzip.open(vcf_path, "rt") as file:
                return self._extract_dimensions_from_opened_vcf(file)
        except (OSError, gzip.BadGzipFile):
            with open(vcf_path, "r") as file:
                return self._extract_dimensions_from_opened_vcf(file)

    def _extract_dimensions_from_opened_vcf(self, file) -> dict:
        """
        Counts the number of samples and variants in the VCF file.
        """
        variant_count = 0
        sample_count = 0
        scaffold_names = set()

        for line in file:
            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
                sample_IDs = header[9:]
                sample_count = len(sample_IDs)
            if not line.startswith("#"):
                variant_count += 1
                scaffold = line.split("\t", 1)[0]
                scaffold_names.add(scaffold)

        return {
            "variants": variant_count,
            "sample_count": sample_count,
            "scaffolds": sorted(list(scaffold_names)),
        }

    def _get_bucket_dimensions_file(self) -> dict:
        """
        If bucket version files exists, download version metadata file from the S3 bucket.
        If doesn't exist, return empty dict.
        """
        try:
            content = self.s3_file_manager.download_s3_file_to_str(
                key=DIMENSIONS_FILE_NAME, bucket_name=self.bucket_name
            )
        except ObjectDoesNotExistError:
            logger.info(f"No bucket versioning file found in the bucket: {self.bucket_name}.")
            return {}
        if not content:
            return {}

        return yaml.safe_load(content)

    def _upload_bucket_dimensions_file(self, version_data: dict) -> None:
        """
        Upload a new version of the .vcf_dimensions.yaml file to the S3 bucket.
        Works for both creating and updating the file.
        """
        text_content = yaml.safe_dump(version_data, sort_keys=False)
        try:
            self.s3_file_manager.upload_str_as_s3_object(
                key=DIMENSIONS_FILE_NAME, content=text_content, bucket_name=self.bucket_name
            )
            logging.info(f"New version updated in the bucket: {self.bucket_name}.")
        except botocore.exceptions.ClientError as e:
            logging.error(f"Failed to upload bucket version file: {e}")

    def download_yaml_to_repo_root(self) -> None:
        """
        Download the .vcf_dimensions.yaml file from the bucket and write it to the repo root.
        """
        # TODO - this method is mostly for dev. consider skipping it. divbase cli download file would achieve the same thing
        yaml_data = self._get_bucket_dimensions_file()
        if not yaml_data:
            print(f"No .vcf_dimensions.yaml found in bucket: {self.bucket_name}.")
            return
        repo_root = Path(__file__).parent.parent.parent
        local_yaml_path = repo_root / ".vcf_dimensions.yaml"
        with open(local_yaml_path, "w") as output_file:
            yaml.safe_dump(yaml_data, output_file, sort_keys=False)
        print(f".vcf_dimensions.yaml downloaded to {local_yaml_path}")

    def get_dimensions_info(self) -> dict:
        """
        Returns the contents of the .vcf_dimensions.yaml file as a dictionary.
        """
        if not self.dimensions_info:
            logger.info("No VCF dimensions index has been created for this bucket as of yet.")
            return {}
        return self.dimensions_info


def create_bucket_manager(project_config: ProjectConfig) -> VCFDimensionIndexManager:
    """
    Helper function to create a BucketVersionManager instance.
    Used by the version and file subcommands of the CLI
    """
    s3_file_manager = create_s3_file_manager(project_config.s3_url)
    return VCFDimensionIndexManager(bucket_name=project_config.bucket_name, s3_file_manager=s3_file_manager)


def show_dimensions_command(project_config: ProjectConfig) -> dict[str, dict]:
    manager = create_bucket_manager(project_config=project_config)
    return manager.get_dimensions_info()
