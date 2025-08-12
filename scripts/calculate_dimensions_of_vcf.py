"""
Helper script that calculates the number of samples and variants in a VCF file and appends the results to a log.
Works with compressed and uncompressed VCF files.

Usage:
    python scripts/calculate_dimensions_of_vcf.py --vcf <path_to_vcf> [--output <path_to_output.tsv>]

    (--output defaults to vcf_dimensions.tsv if not specified)
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
    parser.add_argument(
        "--skip-log",
        action="store_true",
        help="If set, do not write results to the log file. Overridden when --force is also set.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="If set, always calculate and log dimensions, even if the file is already in the log. Overwrites the log entry if it exists. This option ignores --skip-log.",
    )
    return parser.parse_args()


def filename_in_log(filename: str, log_path: Path) -> bool:
    if not log_path.exists():
        return False
    with open(log_path, "r") as f:
        for line in f:
            if line.strip().startswith(filename + "\t"):
                print(f"Previous entry found in log ({log_path}):")
                print("filename\tsample_count\tvariant_count")
                print(f"{line.strip()}\n")
                return True
    return False


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


def write_results_to_log(dimensions: dict, vcf_path: Path, output_path: Path = None, overwrite: bool = False) -> None:
    """
    Write the results to a log file. If overwrite is True, replace any existing entry for this file; this is done by
    going through all existing lines in the log file and carrying all non-matching lines over to a new list; finally,
    the overwritten/updated line is appended to the log.
    """
    header = "filename\tsample_count\tvariant_count\n"
    line = f"{vcf_path.name}\t{dimensions['sample_count']}\t{dimensions['variants']}\n"

    lines = []
    found = False

    if output_path.exists():
        with open(output_path, "r") as f:
            lines = f.readlines()
        if overwrite:
            new_lines = []
            for l in lines:
                if not l.strip().startswith(vcf_path.name + "\t"):
                    new_lines.append(l)
                else:
                    found = True
            lines = new_lines

    if not output_path.exists():
        lines = [header]

    lines.append(line)

    with open(output_path, "w") as f:
        f.writelines(lines)

    if found:
        print(f"Overwrote previous entry for {vcf_path.name} in log at: {output_path}.")
    else:
        print(f"Wrote VCF dimension results to log at: {output_path}.")


def main():
    args = parse_arguments()
    vcf_path = Path(args.vcf)
    output_path = args.output

    skip_log = args.skip_log and not args.force

    if not args.force and filename_in_log(vcf_path.name, output_path):
        print(f"The file {vcf_path.name} is already present in {output_path}. Skipping calculation.")
        print("(run script with --force to rerun calculation and overwrite existing log entry)")
        return

    if vcf_path.name.endswith(".vcf") or vcf_path.name.endswith(".vcf.gz"):
        dimensions = wrapper_calculate_dimensions(vcf_path=vcf_path)
    else:
        print("Invalid file extension. Please provide a .vcf or .vcf.gz file.")
        return

    print(f"Number of samples: {dimensions['sample_count']}")
    print(f"Number of variants: {dimensions['variants']}")

    if not skip_log:
        write_results_to_log(
            dimensions=dimensions,
            vcf_path=vcf_path,
            output_path=output_path,
            overwrite=args.force,
        )
    else:
        print("Skipping writing results to log file, since --skip-log was provided.")


if __name__ == "__main__":
    main()
