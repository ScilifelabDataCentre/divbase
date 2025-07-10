"""
Helper script to generate mock sample metadata from a VCF file.
This script reads a VCF file (either gzipped or uncompressed) and generates a sample metadata file
where each row represents a sample ID, a mock sampling population number, a mock sampling area, a mock sample sex, and the filename that the sample is found in.
"""

import argparse
import gzip
from pathlib import Path
from typing import TextIO


def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate mock sample metadata from one or more VCF files.")
    parser.add_argument(
        "--vcf",
        type=Path,
        required=True,
        help="Comma-separated list of input VCF files (uncompressed or gzipped).",
    )
    return parser.parse_args()


def parse_vcf_file(vcf_path: Path) -> None:
    """
    Function that reads compressed or uncompressed VCF files and passes it onwards to the function that generates the mock sample metadata.
    The code redundancy is intentional to comply with the ruff linter when using context managers for both compressed and uncompressed files.

    """
    print(f"Reading: {vcf_path} ...")
    try:
        with gzip.open(vcf_path, "rt") as file:
            generate_mock_sample_metadata(file, vcf_path)
    except gzip.BadGzipFile:
        with open(vcf_path, "r") as file:
            generate_mock_sample_metadata(file, vcf_path)


def generate_mock_sample_metadata(file: TextIO, vcf_path: Path) -> None:
    """
    Function that generates mock sample metadata from a VCF file. It also counts and displays the number of samples and variants in the VCF file.

    The output file contains the mandatory columns Sample_ID and Filename, as well as three mock columns: Population, Area, and Sex.

    To create some variation across the three mock columns, the three mock columns are generated from lists of different lengths.
    To ensure that the periodicity of the mock area and sex columns are different, the mock sex column dependent on the length of the mock area column.
    Thus, for each sample:
    - The mock area will be "North", "East", "South", "West", and will repeat every 4 samples.
    - The mock population will be a number from 1 to 6, and will repeat every 6 samples.
    - The mock sex will be "F" or "M", and will repeat "F" for 4 samples, then "M" for 4 samples, and so on.
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

    print(f"Number of samples: {sample_count}")
    print(f"Number of variants: {variant_count}")

    mock_area = ["North", "East", "South", "West"]
    mock_population = [1, 2, 3, 4, 5, 6]
    mock_sex = ["F", "M"]

    output_file = f"{vcf_path}_sample_metadata.tsv"

    with open(output_file, "w") as writer:
        writer.write("#Sample_ID\tPopulation\tArea\tSex\tFilename\n")
        for i, sample in enumerate(sample_IDs):
            area = mock_area[i % len(mock_area)]
            population = mock_population[i % len(mock_population)]
            sex = mock_sex[(i // len(mock_area)) % len(mock_sex)]
            writer.write(f"{sample}\t{population}\t{area}\t{sex}\t{vcf_path}\n")
        print(f"Wrote mock sidecar metadata file to: {output_file}")


def main():
    args = parse_arguments()
    vcf_path = args.vcf

    if vcf_path.name.endswith(".vcf") or vcf_path.name.endswith(".vcf.gz"):
        parse_vcf_file(vcf_path)
    else:
        print("Invalid file extension. Please provide a .vcf or .vcf.gz file.")


if __name__ == "__main__":
    main()


# TODO allow for multiple vcf files to be passed in, and generate a single sidecar metadata file for all of them.
