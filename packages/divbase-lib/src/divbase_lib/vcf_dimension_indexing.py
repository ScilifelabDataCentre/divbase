import datetime
import gzip
import logging
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from divbase_lib.exceptions import ObjectDoesNotExistError
from divbase_lib.s3_client import S3FileManager

DIMENSIONS_FILE_NAME = ".vcf_dimensions.yaml"

logger = logging.getLogger(__name__)


@dataclass
class VCFDimensionIndexManager:
    bucket_name: str
    s3_file_manager: S3FileManager

    # not defined till post_init method run.
    dimensions_info: dict = field(init=False)

    def __post_init__(self):
        """
        Loads the .vcf_dimensions.yaml file from the bucket if it exists.
        """
        self.dimensions_info = self._get_bucket_dimensions_file()

    def add_dimension_entry(self, vcf_filename: str) -> None:
        """
        Append a new dimension entry to .vcf_dimensions.yaml if not already present for that VCF file.
        Calls submethods to calculate dimensions and upload the updated YAML file to bucket.

        NOTE! This method does not upload the updated dimensions file to the bucket. The caller must do that
        by calling self._upload_bucket_dimensions_file(self.dimensions_info) after calling this method.
        This way, the caller can add multiple entries and then upload the file only once at the end. Thus,
        the version history of the dimensions file will not be cluttered with many versions that differ only by one entry.
        """

        yaml_data = self.dimensions_info or {"dimensions": []}
        dimensions = yaml_data.get("dimensions", [])

        dimensions = yaml_data.get("dimensions", [])

        vcf_path = Path(vcf_filename)
        vcf_dims = self._wrapper_calculate_dimensions(vcf_path)
        if vcf_dims is None:
            logger.info(f"Skipping indexing of {vcf_filename}: detected DivBase-generated result VCF.")
            return
        timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
        latest_versions_of_bucket_files = self.s3_file_manager.latest_version_of_all_files(bucket_name=self.bucket_name)
        file_version_ID = latest_versions_of_bucket_files.get(vcf_filename, "null")

        existing_entry = next((entry for entry in dimensions if entry.get("filename") == vcf_filename), None)

        if existing_entry and existing_entry.get("file_version_ID_in_bucket") == file_version_ID:
            logger.info(f"Filename '{vcf_filename}' with current version is already present in .vcf_dimensions.yaml.")
            return

        if existing_entry:
            existing_entry["timestamp_from_dimensions_indexing"] = timestamp
            existing_entry["file_version_ID_in_bucket"] = file_version_ID
            existing_entry["dimensions"] = vcf_dims
            logger.info(f"Updated entry for {vcf_filename} in vcf_dimensions index.")
        else:
            new_entry = {
                "filename": vcf_filename,
                "timestamp_from_dimensions_indexing": timestamp,
                "file_version_ID_in_bucket": file_version_ID,
                "dimensions": vcf_dims,
            }
            dimensions.append(new_entry)
            logger.info(f"Added new entry for {vcf_filename} to vcf_dimensions index.")

        yaml_data["dimensions"] = dimensions
        self.dimensions_info = yaml_data

    def remove_dimension_entry(self, vcf_filename: str) -> None:
        """
        Remove a dimension entry from .vcf_dimensions.yaml if present for that VCF file.
        Calls submethods to upload the updated YAML file to bucket.
        """
        yaml_data = self.dimensions_info
        dimensions = yaml_data.get("dimensions", [])

        dimensions = [entry for entry in dimensions if entry.get("filename") != vcf_filename]
        yaml_data["dimensions"] = dimensions

        logger.info(f"Removed entry for {vcf_filename} from .vcf_dimensions.yaml.")

    def get_indexed_filenames(self) -> list[str]:
        """
        Returns a list of all filenames already indexed in .vcf_dimensions.yaml for this bucket.
        """
        yaml_data = self.dimensions_info

        return {
            entry.get("filename"): entry.get("file_version_ID_in_bucket")
            for entry in yaml_data.get("dimensions", [])
            if "filename" in entry
        }

    def _wrapper_calculate_dimensions(self, vcf_path: Path) -> dict:
        """
        Reads compressed or uncompressed VCF files and returns their dimensions.

        Two different entry points to self._extract_dimensions_from_opened_vcf() in
        order to comply with Ruff linting of context managers.
        """
        logger.debug(f"Reading: {vcf_path} ...")
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
            if line.startswith("##DivBase_created"):
                return None
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
            "sample_names": sample_IDs,
        }

    def _get_bucket_dimensions_file(self) -> dict:
        """
        Download VCF dimension metadata file from the S3 bucket.
        Since the file is created in __post_init__, it will exist by the time this method is called.
        But if the upload did not work for whatever reason, there still needs to be error handling for missing file.
        """
        try:
            content = self.s3_file_manager.download_s3_file_to_str(
                key=DIMENSIONS_FILE_NAME, bucket_name=self.bucket_name
            )
        except ObjectDoesNotExistError:
            logger.info(f"No VCF dimensions file found in the bucket: {self.bucket_name}.")
            return {"dimensions": []}
        if not content:
            return {"dimensions": []}

        data = yaml.safe_load(content)
        if not isinstance(data, dict) or "dimensions" not in data:
            logger.warning(
                f"Malformed VCF dimensions file in bucket: {self.bucket_name}. Returning empty dimensions list."
            )
            return {"dimensions": []}
        return data

    def _upload_bucket_dimensions_file(self, dimensions_data: dict) -> None:
        """
        Upload a new version of the .vcf_dimensions.yaml file to the S3 bucket.
        Works for both creating and updating the file.
        """
        text_content = yaml.safe_dump(dimensions_data, sort_keys=False)
        self.s3_file_manager.upload_str_as_s3_object(
            key=DIMENSIONS_FILE_NAME, content=text_content, bucket_name=self.bucket_name
        )
        logging.info(f"New dimensions file uploaded to the bucket: {self.bucket_name}.")

    def get_dimensions_info(self) -> dict:
        """
        Returns the contents of the .vcf_dimensions.yaml file as a dictionary.
        """
        if self.dimensions_info is None or "dimensions" not in self.dimensions_info:
            logger.info("No VCF dimensions have been created for this bucket as of yet.")
            return {"dimensions": []}
        return self.dimensions_info
