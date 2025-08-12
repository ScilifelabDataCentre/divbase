"""
Helper script to generate mock sample metadata from a VCF file to allow for testing the codebase with VCFs that don't have a metadata sidecar file.
This script reads a comma-separated list VCF files (either gzipped or uncompressed) and generates a single sample metadata file based on the sample IDs found in each VCF.
The output file has these columns: sample ID, mock sampling population number, mock sampling area, mock sample sex, and the filename that the sample is found in.
"""

import argparse
import gzip
from pathlib import Path
from typing import TextIO


def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate mock sample metadata from one or more VCF files.")
    parser.add_argument(
        "--vcf",
        type=str,
        required=True,
        help="Comma-separated list of input VCF files (uncompressed or gzipped).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("mock_sample_metadata.tsv"),
        help="Output metadata file (default: mock_sample_metadata.tsv)",
    )
    return parser.parse_args()


def wrapper_get_sample_IDs_from_vcf_file(vcf_path: Path) -> None:
    """
    Function that reads compressed or uncompressed VCF files and passes it onwards to the function that generates the mock sample metadata.
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
    Function that extracts sample IDs from the VCF file header.
    Returns a list of sample IDs.
    """
    for line in file:
        if line.startswith("#CHROM"):
            header = line.strip().split("\t")
            sample_IDs = header[9:]
            return sample_IDs
    return []


def generate_mock_sample_metadata(all_samples: dict[tuple], output_file: Path) -> None:
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

    mock_area = ["North", "East", "South", "West"]
    mock_population = [1, 2, 3, 4, 5, 6]
    mock_sex = ["F", "M"]

    with open(output_file, "w") as writer:
        writer.write("#Sample_ID\tPopulation\tArea\tSex\tFilename\n")
        for i, (sample, vcf_filename) in enumerate(all_samples):
            area = mock_area[i % len(mock_area)]
            population = mock_population[i % len(mock_population)]
            sex = mock_sex[(i // len(mock_area)) % len(mock_sex)]
            writer.write(f"{sample}\t{population}\t{area}\t{sex}\t{vcf_filename}\n")
        print(f"Wrote mock sidecar metadata file to: {output_file}")


def main():
    args = parse_arguments()
    vcf_paths = [Path(vcf_file.strip()) for vcf_file in args.vcf.split(",")]
    output_file = Path(args.output)

    all_samples = []
    for vcf_path in vcf_paths:
        if vcf_path.name.endswith(".vcf") or vcf_path.name.endswith(".vcf.gz"):
            sample_IDs = wrapper_get_sample_IDs_from_vcf_file(vcf_path=vcf_path)
            all_samples.extend((sample, vcf_path.name) for sample in sample_IDs)
        else:
            print("Invalid file extension. Please provide a .vcf or .vcf.gz file.")
    generate_mock_sample_metadata(all_samples=all_samples, output_file=output_file)


if __name__ == "__main__":
    main()
