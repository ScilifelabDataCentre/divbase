from pathlib import Path
import pandas as pd

def tsv_query_command(file: Path, filter: str) -> tuple[pd.DataFrame, str]:
    """
    Query a TSV file for specific key-value pairs. Assumes that the TSV file has a header row 
    and a column named Filename. The user can use any column header as a key, and the values 
    can be any string. Returns the unique filenames that match the query.
    """

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
    query_result = df.query(query_string)

    return query_result, query_string

    #TODO if key not in df.columns, there is currently no warning or error
    #TODO if value not in df[key], there is no warning or error. only if both values are missing 
    # is there is message saying that an empty df is returned
    # TODO what if the user wants to make queries on the Filename column? It should work, but might result in wierd edge-cases?
