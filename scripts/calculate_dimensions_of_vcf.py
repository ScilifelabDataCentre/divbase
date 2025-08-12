"""
Helper script that calculates the number of samples and variants in a VCF file and appends the results to a log.
"""

import argparse
import gzip
from pathlib import Path
from typing import TextIO


def parse_arguments():
    parser = argparse.ArgumentParser(description="Calculate the number of samples and variants in a VCF file.")
    parser.add_argument(
        "--vcf",
        type=str,
        required=True,
        help="Path to the VCF file (uncompressed or gzipped).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("vcf_dimensions.tsv"),
        help="Path to results file with the number of samples and variants (default: vcf_dimensions.tsv)",
    )
    return parser.parse_args()


def wrapper_calculate_dimensions(vcf_path: Path) -> None:
    """
    Wrapper function that reads compressed or uncompressed VCF files and passes it onwards to the function that calculates the dimensions.
    The code redundancy is intentional to comply with the ruff linter when using context managers for both compressed and uncompressed files.

    """
    print(f"Reading: {vcf_path} ...")
    try:
        with gzip.open(vcf_path, "rt") as file:
            return extract_samples_from_opened_vcf(file=file)
    except gzip.BadGzipFile:
        with open(vcf_path, "r") as file:
            return extract_samples_from_opened_vcf(file=file)


def extract_samples_from_opened_vcf(file: TextIO) -> list[str]:
    """
    Function that counts the number of samples (sample columns in the #CHROM header)
    and the number of variants (non-header lines) in the VCF file.
    """
    variant_count = 0
    sample_count = 0

    for line in file:
        if line.startswith("#CHROM"):
            header = line.strip().split("\t")
            sample_IDs = header[9:]
            sample_count = len(sample_IDs)
        if not line.startswith("#"):
            variant_count += 1

    return {
        "variants": variant_count,
        "sample_count": sample_count,
    }


def write_results_to_log(dimensions: dict, vcf_path: Path, output_path: Path = None) -> None:
    """
    Function to write the results to a log file. Checks if the entry already exists.
    """
    header = "filename\tsample_count\tvariant_count\n"
    line = f"{vcf_path.name}\t{dimensions['sample_count']}\t{dimensions['variants']}\n"

    already_logged = False
    if output_path.exists():
        with open(output_path, "r") as f:
            for existing_line in f:
                if existing_line.strip() == line.strip():
                    already_logged = True
                    break

    if not already_logged:
        if not output_path.exists():
            with open(output_path, "w") as f:
                f.write(header)
        with open(output_path, "a") as f:
            f.write(line)
            print(f"Wrote results to log at: {output_path}.")
    else:
        print(f"Entry already exists in log at: {output_path}. Skipping write to log.")


def main():
    args = parse_arguments()
    vcf_path = Path(args.vcf)

    if vcf_path.name.endswith(".vcf") or vcf_path.name.endswith(".vcf.gz"):
        dimensions = wrapper_calculate_dimensions(vcf_path=vcf_path)
    else:
        print("Invalid file extension. Please provide a .vcf or .vcf.gz file.")

    print(f"Number of samples: {dimensions['sample_count']}")
    print(f"Number of variants: {dimensions['variants']}")

    write_results_to_log(dimensions=dimensions, vcf_path=vcf_path, output_path=args.output)


if __name__ == "__main__":
    main()
