"""
Helper script to generate mock sample metadata from a VCF file to allow for testing the codebase with VCFs that don't have a metadata sidecar file.
This script reads a comma-separated list VCF files (either gzipped or uncompressed) and generates a single sample metadata file based on the sample IDs found in each VCF.
The output file has these columns: sample ID, mock sampling population number, mock sampling area, mock sample sex, and the filename that the sample is found in.

Usage:
python scripts/generate_mock_sample_metadata.py --vcf {vcf_filename} --output {metadata_filename} --columns {number_of_random_columns_to_add} --add-warning-column
"""

import argparse
import gzip
import random
import string
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
    parser.add_argument(
        "--columns",
        type=int,
        default=3,
        help="Number of additional random columns to generate after the first 3 legacy columns (Population, Area, Sex). Default: 3",
    )
    parser.add_argument(
        "--add-warning-column",
        action="store_true",
        help="Add a column that will trigger a warning in the TSV validator (e.g., array notation)",
    )
    parser.add_argument(
        "--add-error-column",
        action="store_true",
        help="Add a column that will trigger a hard error in the TSV validator (mixed-type array)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=12345,
        help="Random seed for reproducibility (default: 12345)",
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


def generate_random_values(n):
    vals = []
    for _ in range(n):
        t = random.choice(["int", "float", "str"])
        if t == "int":
            vals.append(str(random.randint(1, 100)))
        elif t == "float":
            vals.append(f"{random.uniform(1, 100):.2f}")
        else:
            vals.append(random_string(5))
    return vals


def generate_mock_sample_metadata(
    all_samples: list[str], output_file: Path, num_columns: int, seed: int, add_warning: bool, add_error: bool
) -> None:
    """
    Generate mock sample metadata with legacy columns (Population, Area, Sex), user-specified number of random columns, and optional warning column.

    The "legacy" columns were used in the first iteration of this script, and are kept for backwards compatibility.
    """
    random.seed(seed)
    legacy_cols = ["Population", "Area", "Sex"]
    random_col_names = [f"Col{i + 1}" for i in range(num_columns)]
    col_names = legacy_cols[:]
    if add_error:
        col_names.append("ErrorCol")
    if add_warning:
        col_names.append("WarningCol")
    col_names += random_col_names

    mock_area = ["North", "East", "South", "West"]
    mock_population = [1, 2, 3, 4, 5, 6]
    mock_sex = ["F", "M"]

    with open(output_file, "w") as writer:
        writer.write("#Sample_ID\t" + "\t".join(col_names) + "\n")
        for i, sample in enumerate(all_samples):
            # Legacy columns
            area = mock_area[i % len(mock_area)]
            population = mock_population[i % len(mock_population)]
            sex = mock_sex[(i // len(mock_area)) % len(mock_sex)]
            row = [str(population), area, sex]
            if add_error:
                # Always generate a mixed-type value (random order, using the seed)
                parts = []
                for _ in range(3):
                    if random.choice([True, False]):
                        parts.append(str(random.randint(1, 100)))
                    else:
                        parts.append(f'"{random_string(random.randint(4, 6))}"')
                error_col_val = f"[{', '.join(parts)}]"
                row.append(error_col_val)
            if add_warning:
                # Always generate a bracketed list of strings (not mixed-type)
                parts = [f'"{random_string(random.randint(4, 6))}"' for _ in range(3)]
                warning_col_val = f"[{', '.join(parts)}]"
                row.append(warning_col_val)
            row += generate_random_values(num_columns)
            writer.write(f"{sample}\t" + "\t".join(row) + "\n")
        print(f"Wrote mock sidecar metadata file to: {output_file}")


def random_string(length=6):
    return "".join(random.choices(string.ascii_letters, k=length))


def main():
    args = parse_arguments()
    vcf_paths = [Path(vcf_file.strip()) for vcf_file in args.vcf.split(",")]
    output_file = Path(args.output)
    num_columns = args.columns
    seed = args.seed
    add_warning = args.add_warning_column
    add_error = args.add_error_column

    all_samples = []
    for vcf_path in vcf_paths:
        if vcf_path.name.endswith(".vcf") or vcf_path.name.endswith(".vcf.gz"):
            sample_IDs = wrapper_get_sample_IDs_from_vcf_file(vcf_path=vcf_path)
            all_samples.extend(sample for sample in sample_IDs)
        else:
            print("Invalid file extension. Please provide a .vcf or .vcf.gz file.")
    generate_mock_sample_metadata(
        all_samples=all_samples,
        output_file=output_file,
        num_columns=num_columns,
        seed=seed,
        add_warning=add_warning,
        add_error=add_error,
    )


if __name__ == "__main__":
    main()
