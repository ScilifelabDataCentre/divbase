import gzip
import logging
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class VCFDimensionCalculator:
    """
    Calculates dimensions (samples, variants, scaffolds) from VCF files.

    This is a pure utility class with no side effects - it only reads VCF files
    and returns their dimensions. Database operations are handled via the API.
    """

    def calculate_dimensions(self, vcf_path: Path) -> dict | None:
        """
        Calculate dimensions from a VCF file.

        Args:
            vcf_path: Path to VCF file (can be .vcf or .vcf.gz)

        Returns:
            Dict with keys: variants, sample_count, scaffolds, sample_names
            None if this is a DivBase-generated result file (should be skipped)
        """
        logger.debug(f"Reading: {vcf_path} ...")
        try:
            with gzip.open(vcf_path, "rt") as file:
                return self._extract_dimensions_from_opened_vcf(file)
        except (OSError, gzip.BadGzipFile):
            with open(vcf_path, "r") as file:
                return self._extract_dimensions_from_opened_vcf(file)

    def _extract_dimensions_from_opened_vcf(self, file) -> dict | None:
        """
        Parse VCF file and extract dimensions.

        Returns None if this is a DivBase-generated result file.
        """
        variant_count = 0
        sample_count = 0
        scaffold_names = set()
        sample_IDs = []

        for line in file:
            # Skip DivBase-generated result files
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


# class VCFDimensionIndexManager:
#     """Class that handles API calls to VCF dimensions database table."""

#     def __init__(self, bucket_name: str, api_base_url: str = None, auth_token: str = None):
#         self.bucket_name = bucket_name
#         self.api_base_url = os.environ.get(
#             "DIVBASE_API_URL", "http://fastapi:8000/api"
#         )  # TODO there might be a better way of passing this. this works with the tests but will run once for every instance
#         self.auth_token = auth_token
#         self._project_id = None  # init to None to avoid API call in init method

#     def _ensure_project_id(self) -> int:
#         """
#         Fetch project_id from API if not already fetched.
#         """
#         if self._project_id is not None:
#             return self._project_id

#         response = httpx.get(
#             f"{self.api_base_url}/v1/vcf-dimensions/lookup/project-by-bucket/{self.bucket_name}",
#             headers=self._auth_headers(),
#         )

#         if response.status_code == 404:
#             raise ValueError(f"Project not found for bucket: {self.bucket_name}")

#         response.raise_for_status()
#         self._project_id = response.json()["project_id"]
#         return self._project_id

#     # def create_or_update_vcf_metadata(self, db: Session, vcf_metadata: dict) -> None:
#     #     """
#     #     Create or update VCF metadata entry via API.

#     #     Note: vcf_metadata dict must include "project_id" key.
#     #     """
#     #     entry = create_or_update_vcf_metadata(db, vcf_metadata)

#     #     logger.info(f"VCF metadata created/updated for {entry.vcf_file_s3_key} in project {entry.project_id}")

#     # def create_or_update_skipped_vcf(self, skipped_data: dict) -> dict:
#     #     """Mark a VCF as skipped via API."""
#     #     response = httpx.post(
#     #         f"{self.api_base_url}/v1/vcf-dimensions/create-skipped",
#     #         json=skipped_data,
#     #         headers=self._auth_headers(),
#     #     )
#     #     response.raise_for_status()
#     #     return response.json()

#     # def get_dimensions_info(self) -> dict:
#     #     """Get all dimensions for this project."""
#     #     project_id = self._ensure_project_id()

#     #     response = httpx.get(
#     #         f"{self.api_base_url}/v1/vcf-dimensions/list/project/{project_id}",
#     #         headers=self._auth_headers(),
#     #     )
#     #     if response.status_code == 404:
#     #         return {"vcf_files": []}
#     #     response.raise_for_status()

#     #     return self._format_for_compatibility(response.json())

#     # def get_skipped_files(self) -> dict:
#     #     """Get all skipped VCF files for this project."""
#     #     try:
#     #         project_id = self._ensure_project_id()
#     #     except ValueError:
#     #         return {}

#     #     with SyncSessionLocal() as db:
#     #         entries = get_skipped_vcfs_by_project_worker(db, project_id)

#     # skipped_data = {
#     #     "project_id": project_id,
#     #     "skipped_file_count": len(entries),
#     #     "skipped_files": [
#     #         {
#     #             "vcf_file_s3_key": entry.vcf_file_s3_key,
#     #             "s3_version_id": entry.s3_version_id,
#     #         }
#     #         for entry in entries
#     #     ],
#     # }

#     # return {entry["vcf_file_s3_key"]: entry["s3_version_id"] for entry in skipped_data.get("skipped_files", [])}

#     # def delete_vcf_metadata(self, vcf_file_s3_key: str, project_id: int) -> None:
#     #     """Delete VCF metadata entry."""
#     #     response = httpx.delete(
#     #         f"{self.api_base_url}/v1/vcf-dimensions/delete/project/{project_id}/file/{vcf_file_s3_key}",
#     #         headers=self._auth_headers(),
#     #     )
#     #     response.raise_for_status()

#     # def delete_skipped_vcf(self, vcf_file_s3_key: str, project_id: int) -> None:
#     #     """Delete skipped VCF entry."""
#     #     response = httpx.delete(
#     #         f"{self.api_base_url}/v1/vcf-dimensions/delete-skipped/project/{project_id}/file/{vcf_file_s3_key}",
#     #         headers=self._auth_headers(),
#     #     )
#     #     response.raise_for_status()

#     def _auth_headers(self) -> dict:
#         if self.auth_token:
#             return {"Authorization": f"Bearer {self.auth_token}"}
#         return {}

#     def _format_for_compatibility(self, api_response: dict) -> dict:
#         """Convert API response to old format for backward compatibility."""

#         # TODO it would make more sense to keep the same keys throughout. but the tests historically expect these keys, so reformat for now
#         vcf_files_list = []
#         for entry in api_response.get("vcf_files", []):
#             vcf_files_list.append(
#                 {
#                     "vcf_file_s3_key": entry["vcf_file_s3_key"],
#                     "s3_version_id": entry["s3_version_id"],
#                     "sample_names": entry.get("samples", []),
#                     "scaffolds": entry.get("scaffolds", []),
#                     "variant_count": entry.get("variant_count", 0),
#                     "sample_count": entry.get("sample_count", 0),
#                     "indexed_at": entry.get("indexed_at"),
#                 }
#             )
#         return {"vcf_files": vcf_files_list}


# def create_vcf_dimension_manager(bucket_name: str, auth_token: str = None) -> VCFDimensionIndexManager:
#     """Factory function to create VCFDimensionIndexManager."""
#     return VCFDimensionIndexManager(bucket_name=bucket_name, auth_token=auth_token)
