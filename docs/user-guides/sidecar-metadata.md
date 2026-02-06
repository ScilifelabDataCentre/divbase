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

#### Mandatory content

1. The first row must be a header row and the first column must be named `Sample_ID`.
2. The `Sample_ID` column must contain the exact names of the samples as they are spelled in the VCF files. This will already be handled if user has run a `divbase-cli dimensions update` job and, after its completion, has generated a pre-filled template with: `divbase-cli dimensions create-metadata-template`
3. The `Sample_ID` column can only contain one sample name per row. This is different from the user-defined columns that can take arrays of values for each cell in a column.
4. Every column need to be tab separated for all rows.

#### User-defined columns

After the `Sample_ID` column has been populated, users can add any columns and values to the TSV.

1. It is the user's responsibility to ensure that the spelling of column headers and values is consistent. When filtering on the sidecar metadata, the exact spelling must be used for the filters.
2. The user-defined columns can be either numeric or string type. Try to avoid mixing string and numeric values in the same column is possible. If a mix of string and numerical data is used in the same column, the system will treat them all as strings, which might lead to unexpected filtering results when running queries. The DivBase backend uses [`Pandas`](https://pandas.pydata.org/) to automatically infer column type based on its data, so there is no need to specify in the TSV whether the values is numerical or string.
3. Use English decimal notation (.) and not comma (,) when entering decimals. This ensures that the data is correctly loaded by `Pandas`.
4. Semicolon-separated values are supported in TSV cells to represent arrays of values. This allows users to have samples that can belong to multiple values in the same column. For instance belong to two different groups or categories. This works with both numerical and string data.

#### Example

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

### Overview: querys are applied as filters on columns in the TSV

Queries on the sidecar sample metadata TSV can be done with the `divbase-cli query tsv` command. The filters that the user want to query on needs entered as a string (i.e. enclosed in quotes, `""`).

The TSV query syntax is `"Key1:Value1,Value2;Key2:Value3,Value4"`, where `Key1:`...`Key2:` are the column header names in the TSV, and `Value1`...`Value4` are the values. Multiple filter values for a key are separated by commas, and multiple keys are separated by semicolons. There can be any number keys and values to filter on, but it is up to the user to write queries that return useful results.

!!! note
    Please note that semicolon (`;`) is used for different purposes in the TSV (multi-value cells) and in the query syntax (perform queries on multiple columns)!

Filtering is inclusive by default. This applies both for the filter values and the cell values:

- If a filter contains multiple values, e.g. `"Area:North,West"`, the row is included if at least one of the filter values matches any value in the cell. I.e. a row with `North`, and a row with `West` will both be returned from this filter.
- If a cell in the TSV contains multiple values separated by a semicolon as explained in [User-defined columns](#user-defined-columns) (e.g., `North;West`), the row is included if any of those values match the filter. Filters with `"Area:North"`, `"Area:West"`, and `"Area:North,West"` will all return the row with the array value `North;West`.

For example, if the user wants to query the TSV on column `Area` for all samples that contain the value `North`,:

```bash
divbase-cli query tsv "Area:North"
```

It is also possible to run a sidecar sample metadata query as part of a VCF query by adding the query as a sting to the flag `--tsv-filter`:

```bash
divbase-cli query bcftools-pipe --tsv-filter "Area:North" --command "view -s SAMPLES"
```

Please also see the documentation on [DivBase Query Syntax for VCF data](query-syntax.md) for more details on how that command works.

!!! note
    To reiterate what was written in the [User-defined columns](#user-defined-columns) section above: it the user's responsibility to ensure that the spelling of column headers and values is consistent. When filtering on the sidecar metadata, the exact spelling must be used for the filters.

### Filtering on string columns

Queries on string columns are straight-forward in the sense that each semicolon-separated value in the TSV are treated as discrete values.

As explained above, commas can be used to write multi-values filters. For instance, the query:

```bash
divbase-cli query tsv "Area:North,South,East"
```

will return all samples where **at least one** of the semicolon-separated values in the Area column matches any of the filter values (`North`, `South`, or `East`).

Comma-separated: "Area:North,South,East"
OR logic: if ANY cell value matches ANY filter value, the row matches
Example: "key2:value3,value4" matches cells containing value3, value4, value3;value4, or value4;value3

### Filtering on numerical columns

For numerical columns, it is possible to filter on the following operations:

- Inequalities:
  - Examples: `"Weight:>25"` or `"Weight:>=20,<=40"` or `"Weight:<100"`.
  - Note" The inequality operator must be expressed relative to the key, i.e. `"Weight:>25"`. The reverse notation `"Weight:25<"` is not supported.
  - The syntax only accepts `<=` and `>=` since this is the syntax of Python. The forms =< and => are not accepted and will return an error.
- Range (inclusive):
  - Example: `"Weight:20-40"`
- Discrete values:
  - Example: `"Weight:25,30,35"`

Furthermore, it is possible to combine filters on inequalities, ranges, and discrete values using inclusive OR logic. This means that if any one of the specified conditions is satisfied for a cell, the row will be included in the results. For example:

- `"Weight:<2,4"` returns rows where the value is less than 2 **or** equal to 4
- `"Weight:1-2,4"` returns rows where the value is in the range 1–2 **or** equal to 4
- `"Weight:>5,1-2,4"` returns rows where the value is greater than 5 **or** in the range 1–2 **or** equal to 4
- `"Weight:>10,<2,5-7"` returns rows where the value is greater than 10 **or** less than 2 **or** in the range 5–7

## Trying out a query

```bash
divbase-cli query tsv "Area:Northern Portugal"
```

- [TO BE IMPLEMENTED] what to do if a query references a column that does not exist. E.g. `divbase-cli query tsv "Area:Northern Portugal"` when Area does not exist? This should probably give a warning and not just return nothing

- [TO BE IMPLEMENTED] what to do if a query references a column value. E.g. `divbase-cli query tsv "Area:Northern Portugal"` when Northern Portugal does not exist in the column? This should probably also give a warning and not just return nothing, but nothing is a result here and not a syntax problem...
