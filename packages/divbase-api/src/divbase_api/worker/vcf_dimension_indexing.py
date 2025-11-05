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
