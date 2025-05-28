"""
Query subcommand for the divbase_tools CLI.
"""

from pathlib import Path
import typer
from rich import print
import pandas as pd


query_app = typer.Typer(
    help="Query the metadata for the VCF files stored in the bucket.", no_args_is_help=True
)

@query_app.command("tsv")
def tsv_query_command(
    file: Path = typer.Option(
        default=Path("./sample_metadata.tsv"),
        help="Path to the tsv metadata file.",
    ),
    filter: str = typer.Option(None, 
        help="""String consisting of keys:values in the tsv file to filter on.
        The syntax is 'Key1:Value1,Value2;Key2:Value3,Value4', where the key
        are the column header names in the tsv, and values are the column values. 
        Multiple values for a key are separated by commas, and multiple keys are 
        separated by semicolons.
        E.g. 'Area:West of Ireland,Northern Portugal;Sex:F'.""",
        ),
    show_sample_results: bool = typer.Option(
        default=False,
        help="Print sample_ID and Filename results from the query.",
    ),
):
    """
    Query a TSV file for specific key-value pairs. Assumes that the TSV file has a header row 
    and a column named Filename. The user can use any column header as a key, and the values 
    can be any string. Returns the unique filenames that match the query.
    """
    print(f"Parsing: {file}\n")

    df = pd.read_csv(file, sep="\t")
    df.columns = df.columns.str.lstrip("#")

    key_values = filter.split(";")
    filter_conditions = []
    for key_value in key_values:
        key, values = key_value.split(":", 1)
        values_list = values.split(",")

        if key in df.columns:
            condition = f"{key} in {values_list}"
            filter_conditions.append(condition)

    query_string = " and ".join(filter_conditions)
    result = df.query(query_string)
    
    if show_sample_results:
        print(f"Name and file for each sample in query results:")  
        print(f"{result[["Sample_ID", "Filename"]].to_string(index=False)}\n")

    unique_filenames = result["Filename"].unique()
    print(f"Unique filenames for query {query_string}: {unique_filenames}")    

    #TODO if key not in df.columns, there is currently no warning or error
    #TODO if value not in df[key], there is no warning or error. only if both values are missing 
    # is there is message saying that an empty df is returned
    # TODO what if the user wants to make queries on the Filename column? It should work, but might result in wierd edge-cases?
