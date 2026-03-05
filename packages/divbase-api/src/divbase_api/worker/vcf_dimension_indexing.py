import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List

logger = logging.getLogger(__name__)


@dataclass
class VCFDimensions:
    variants: int
    sample_count: int
    scaffolds: List[str]
    sample_names: List[str]


class VCFDimensionCalculator:
    """
    Calculates dimensions (samples, variants, scaffolds) from VCF files.

    This is a pure utility class with no side effects - it only reads VCF files
    and returns their dimensions. Database operations are handled via the API.
    """

    def calculate_dimensions(self, vcf_path: Path) -> VCFDimensions | None:
        """
        Calculate dimensions from a VCF file using bcftools and CSI indexing.

        Returns None if this is a DivBase-generated result file (should be skipped).
        This is checked for when the header is parsed in _extract_sample_names_from_vcf_header.

        Note: bcftools index --csi requires bgzipped input. Plain .vcf files are bgzipped to a
        temporary .vcf.gz file for indexing, then the temporary file is cleaned up internally.
        """
        logger.debug(f"Reading: {vcf_path} ...")

        sample_names = self._extract_sample_names_from_vcf_header(vcf_path)
        if sample_names is None:
            return None

        indexing_path = vcf_path

        if vcf_path.suffix == ".vcf":
            bgzipped_temp = Path(str(vcf_path) + ".gz")
            self._bgzip_vcf(vcf_path, bgzipped_temp)
            indexing_path = bgzipped_temp

        try:
            csi_index_path = self._index_vcf_with_csi(vcf_path=indexing_path)
            scaffold_names, variant_count = self._extract_scaffold_names_and_variant_count_from_csi_index(
                csi_index_path=csi_index_path
            )
        finally:
            # Clean up temp bgzipped file and its CSI index; these are not tracked by tasks.py
            if bgzipped_temp is not None:
                bgzipped_temp_csi = Path(str(bgzipped_temp) + ".csi")
                if bgzipped_temp.exists():
                    bgzipped_temp.unlink()
                if bgzipped_temp_csi.exists():
                    bgzipped_temp_csi.unlink()

        return VCFDimensions(
            variants=variant_count,
            sample_count=len(sample_names),
            scaffolds=sorted(scaffold_names),
            sample_names=sample_names,
        )

    def _extract_sample_names_from_vcf_header(self, vcf_path: Path) -> List[str] | None:
        """
        Extract sample names from the VCF header using bcftools. Reads only the header lines, so this is fast even for large VCFs.

        Returns None if the VCF turns out to be a is a DivBase-generated result file (Custom header "##DivBase_created").
        """
        try:
            result = subprocess.run(
                ["bcftools", "view", "-h", str(vcf_path)],
                check=True,
                stdout=subprocess.PIPE,
                text=True,
                stderr=subprocess.DEVNULL,  # Important! Would otherwise need to parse out potential stderr downstream
            )
            for line in result.stdout.splitlines():
                if line.startswith("##DivBase_created"):
                    return None
                if line.startswith("#CHROM"):
                    cols = line.strip().split("\t")
                    sample_names = cols[9:] if len(cols) > 9 else []
                    logger.info(f"Extracted {len(sample_names)} sample names from the VCF header.")
                    return sample_names
            return []
        except subprocess.CalledProcessError as e:
            logger.error(f"Error extracting sample names from the VCF header {vcf_path}: {e}")
            raise

    def _bgzip_vcf(self, vcf_path: Path, output_path: Path) -> None:
        """
        Bgzip a plain .vcf file to a bgzipped .vcf.gz file using bcftools view.
        Required because bcftools index --csi only accepts bgzipped input.
        """
        try:
            subprocess.run(
                ["bcftools", "view", "-Oz", "-o", str(output_path), str(vcf_path)],
                check=True,
                stderr=subprocess.DEVNULL,
            )
            logger.info(f"Bgzipped {vcf_path} to {output_path} for CSI indexing.")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error bgzipping {vcf_path}: {e}")
            raise

    def _index_vcf_with_csi(self, vcf_path: Path) -> Path:
        """
        Index the VCF file with CSI index using bcftools.
        The CSI index is then used to extract dimensions information.
        Input must be bgzipped (.vcf.gz); use _bgzip_vcf first for plain .vcf files.
        """
        csi_index_path = vcf_path.with_suffix(vcf_path.suffix + ".csi")
        try:
            subprocess.run(
                ["bcftools", "index", "--csi", str(vcf_path), "-o", str(csi_index_path)],
                check=True,
                stderr=subprocess.DEVNULL,
            )
            logger.info(f"Successfully indexed {vcf_path} with a CSI index.")
            return csi_index_path
        except subprocess.CalledProcessError as e:
            logger.error(f"Error indexing {vcf_path} with CSI index: {e}")
            raise

    def _extract_scaffold_names_and_variant_count_from_csi_index(self, csi_index_path: Path) -> tuple[list[str], int]:
        """
        Extract scaffold names and total variant count from the CSI index file using bcftools.

        The output of bcftools index --stats have three columns: contig name, contig length (. if unknown) and number of records for the contig
        """
        try:
            result = subprocess.run(
                ["bcftools", "index", "--stats", str(csi_index_path)],
                check=True,
                stdout=subprocess.PIPE,
                text=True,
                stderr=subprocess.DEVNULL,
            )
            scaffold_names = []
            variant_count = 0
            for line in result.stdout.strip().split("\n"):
                cols = line.split("\t")
                if len(cols) >= 3:
                    scaffold_names.append(cols[0])
                    variant_count += int(cols[2])
            logger.info(
                f"Extracted {len(scaffold_names)} scaffold names and {variant_count} variants from the CSI index."
            )
            return scaffold_names, variant_count
        except subprocess.CalledProcessError as e:
            logger.error(f"Error extracting scaffold names/variant count from the CSI index {csi_index_path}: {e}")
            raise
