"""
Script to load a sample metadata TSV file to a Pandas DataFrame using the SharedMetadataValidator and
display the df in the terminal. Can be used to ensure that the TSV results in the expected DataFrame structure,
and to inspect the DataFrame for debugging purposes.

Usage:

python scripts/tsv_to_dataframe.py --tsv path/to/metadata.tsv

"""

import argparse
from pathlib import Path

import pandas as pd

from divbase_lib.metadata_validator import SharedMetadataValidator


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Load a sample metadata TSV file to Pandas dataframe using the SharedMetadataValidator."
    )
    parser.add_argument(
        "--tsv",
        type=str,
        required=True,
        help="Path to the TSV file to be loaded.",
    )
    return parser.parse_args()


def tsv_to_dataframe(file_path) -> pd.DataFrame | None:
    """
    Reads a TSV file and returns a pandas DataFrame. Just runs the loading and validation logic, but does not
    print the results like the client-side ClientSideMetadataTSVValidator does.

    Allows for inspection of of the dataframe.
    """
    validator = SharedMetadataValidator(file_path=Path(file_path), project_samples=set(), skip_dimensions_check=True)
    result = validator.load_and_validate()
    return result.df


def main():
    args = parse_arguments()
    df = tsv_to_dataframe(args.tsv)
    if df is not None:
        print(df.head())
    else:
        print("Failed to load DataFrame. Check if the file exists and is a valid TSV.")

    pd.set_option("display.max_columns", None)  # Show all columns
    pd.set_option("display.width", 120)  # Set terminal width
    pd.set_option("display.max_rows", None)  # Show all rows

    print(df)


if __name__ == "__main__":
    main()
