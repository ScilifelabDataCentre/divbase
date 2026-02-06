# Sidecar Metadata TSV files: creating and querying sample metadata files

TODO

- rationale: what is this and how can it be used

## Creating a sidecar TSV for a DivBase project

If the dimensions VCF files in the project have been cached in DivBase, a template metadata file with the sample names pre-filled can be created with:

```bash
divbase-cli dimensions create-metadata-template
```

Note! there can be multiple TSVs in the same project and it is possible to call them for the queries with the `--metadata-tsv-name` flag.

### Sidecar TSV format requirements

**Mandatory content:**

1. The first row must be a header row and the first column must be named `Sample_ID`.
2. The `Sample_ID` column must contain the exact names of the samples as they are spelled in the VCF files. This will already be handled if user has run a `divbase-cli dimensions update` job and, after its completion, has generated a pre-filled template with: `divbase-cli dimensions create-metadata-template`
3. The `Sample_ID` column can only contain one sample name per row. This is different from the user-defined columns that can take arrays of values for each cell in a column.
4. Every column need to be tab separated for all rows.

**User-defined columns:**

After the `Sample_ID` column has been populated, users can add any columns and values to the TSV.

1. It is the user's responsibility to ensure that the spelling of column headers and values is consistent. When filtering on the sidecar metadata, the exact spelling must be used for the filters.
2. The user-defined columns can be either numeric or string type. Try to avoid mixing string and numeric values in the same column is possible. If a mix of string and numerical data is used in the same column, the system will treat them all as strings, which might lead to unexpected filtering results when running queries. The DivBase backend uses [`Pandas`](https://pandas.pydata.org/) to automatically infer column type based on its data, so there is no need to specify in the TSV whether the values is numerical or string.
3. Use English decimal notation (.) and not comma (,) when entering decimals. This ensures that the data is correctly loaded by `Pandas`.
4. Semicolon-separated values are supported in TSV cells to represent arrays of values. This allows users to have samples that can belong to multiple values in the same column. For instance belong to two different groups or categories. This works with both numerical and string data.

**Example:**

This example illustrates how a sidecar sample metadata TSV can look like. The mandatory requirement are fulfilled (heading, `Sample_ID` column, tab-separated file). The user-defined column contains examples of a numerical column (`Population`) and a string column (`Area`). In some cells, semicolons (`;`) are used to assign multiple values to the same sample and column.

```text
#Sample_ID Population Area
S1 1 North
S2 2;4 East
S3 3 West;South
S4 4 West
S5 5 North
S6 6 East
S7 1;3;5 South
S8 2 West
```

TODOs:

- [TO BE IMPLEMENTED] consider changing the mandatory column name from `Sample_ID` to `Sample`
- [TO BE IMPLEMENTED] what happens if a TSV does not contain all the samples in the DivBase project? There should probably be a warning, but not an error?
- [TO BE IMPLEMENTED] what happens if a sample name is misspelled in the TSV? a warning? can this be checked against the dimensions show?
- [TO BE IMPLEMENTED] what happens if a sample is duplicated in the file. what happens if the sample name is duplicated but not the values (diverging duplicate)?

## Query Syntax for sidecar metadata

from `divbase-cli query tsv -h` docstring:
String consisting of keys:values in the tsv file to filter on. The syntax is 'Key1:Value1,Value2;Key2:Value3,Value4', where the key are the column header names in the tsv, and values are the column values. Multiple values for a key are separated by commas, and multiple keys are separated by semicolons. When multple keys are provided, an intersect query will be performed. E.g. 'Area:West of Ireland,Northern Portugal;Sex:F'.

- [TO BE IMPLEMENTED] filtering based on ranges (ints and maybe floats) and not just on strings, e.g. in range 31 - 50. etc...
- [TO BE IMPLEMENTED] add more Set Operations: union, intersection, difference, symmetric difference. Be clear on the default behaviour

Please do not mix numerical and string values in the same column!

For numeric columns, you can filter on:

- Inequalities: 'Weight:>25' or "Weight:>=20,<=40" or "Weight:<100". The inequality operator must be expressed relative to the Key, i.e. for 'Weight:>25' the reverse notation 'Weight:25<' is not supported.
- Range (inclusive): 'Weight:20-40'
- Discrete values: 'Weight:25,30,35'

The syntax only accepts `<=` and `>=` since this is the syntax of Python. The forms =< and => are not accepted and will return an error.

It is possible to combine filters on inequalities, ranges, and discrete values to an OR logic if desired. For example:

Weight:<2,4 → values less than 2 OR equal to 4
Weight:1-2,4 → values in range 1-2 OR equal to 4
Weight:>5,1-2,4 → values greater than 5 OR in range 1-2 OR equal to 4
Weight:>10,<2,5-7 → values >10 OR <2 OR in range 5-7

note that semicolon is allowed in in cells in the TSV, but have another meaning in the query syntax!

TODO write pytests that ensure that these numerical filters work

## Trying out a query

```bash
divbase-cli query tsv "Area:Northern Portugal"
```

- [TO BE IMPLEMENTED] what to do if a query references a column that does not exist. E.g. `divbase-cli query tsv "Area:Northern Portugal"` when Area does not exist? This should probably give a warning and not just return nothing

- [TO BE IMPLEMENTED] what to do if a query references a column value. E.g. `divbase-cli query tsv "Area:Northern Portugal"` when Northern Portugal does not exist in the column? This should probably also give a warning and not just return nothing, but nothing is a result here and not a syntax problem...
