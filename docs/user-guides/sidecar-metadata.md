# Sidecar Metadata TSV files: creating and querying sample metadata files

TODO

- rationale: what is this and how can it be used

assumes that update dimensions has been run for the latest data

!!! Notes
    There is a CLI command to help check that a user-defined sample metadata TSV file aligns with the requirements described on this page. This validator tool will be [described in its own section below](#validating-a-sidecar-metadata-tsv-with-divbase-cli), but, in short, it can be run with:

    ```bash
    divbase-cli dimensions validate-metadata-file path/to/your/sample_metadata.tsv
    ```

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

TODO non empty, unique, no duplication, no semicolones in Sample_ID

#### User-defined columns

After the `Sample_ID` column has been populated, users can add any columns and values to the TSV.

!!! Warning
    It is the user's responsibility to ensure that the spelling of column headers and values is consistent. When filtering on the sidecar metadata, the exact spelling must be used for the filters. This includes matching upper and lower case letters.

To ensure that user-defined metadata can be used in DivBase, we ask you follow the following constraints and considerations:

1. The user-defined columns can be either numeric or string type. Mixing string and numeric values in the same column is not allowed; if a mix is detected, DivBase will raise an error and reject the file. The DivBase backend uses [`Pandas`](https://pandas.pydata.org/) to automatically infer column type based on its data, so there is no need to specify in the TSV whether the values are numerical or string.
2. Commas are not supported for the TSV and the DivBase system will send an error message if it detects any TSV cells with commas in them. Commas can have different meanings in different notation systems and to avoid confusion and to keep it simple, DivBase will simply not handle commas. Note that commas are used in the [Query syntax](#query-syntax-for-sidecar-metadata) for a different purpose. For decimals, use English decimal notation (.) and not comma (,). DivBase allows one single delimiter for enumerations in the TSV files and that is the semicolon (;) as will be described in the bullet.
3. Semicolon-separated values are supported in TSV cells to represent arrays of values. This allows users to have samples that can belong to multiple values in the same column. For instance belong to two different groups or categories. This works with both numerical and string data (e.g. "2;4;21" or "North; North-West"). Note that this might make the process of writing queries on the more complex than if just a single value is use for each cell.
4. As outlined above, the only characters with special meaning or restrictions in the TSV are `#`, `,`, `;`, and `\t` (tab). Other special characters should be supported, but please be aware that Your Milage May Vary. Some common cases that have been tested and are supported include hyphens (`-`), e.g.`North-West`),   diacritic unicodecharacters like `å`,`ä`,`ö`.
5. Leading and trailing whitespaces are removed by the DivBase backend in order to ensure robust filtering and pattern matching. Whitespaces inside strings will be preserved. For instance: " Sample 1 " will be processed as "Sample 1".

TODO - add info on No duplicate column names, no empty column names

#### Example

This example illustrates how a sidecar sample metadata TSV can look like. The mandatory requirement are fulfilled (heading, `Sample_ID` column, tab-separated file). The user-defined column contains examples of a numerical column (`Population`) and a string column (`Area`). In some cells, semicolons (`;`) are used to assign multiple values to the same sample and column.

```text
#Sample_ID Population Area Weight
S1 1 North 12.1
S2 2;4 East 18.8
S3 3 West;South 15.0
S4 4 West 20.2
S5 5 North 16.1
S6 6 East 25.2
S7 1;3;5 South 22.6
S8 2 West 19.5
```

### Validating a sidecar metadata TSV with `divbase-cli`

Manually checking that a TSV fulfills the DivBase requirement can be tedious. To help users validate their sidecar TSV files, the following CLI command has been implemented:

```bash
divbase-cli dimensions validate-metadata-file path/to/your/sample_metadata.tsv
```

The validation runs on the users local computer and not as a job on the DivBase server. It is intendend to be used on sidecar metadata TSV files before they are uploaded to the DivBase project. The validator will check the formatting requirements as described in [Mandatory contents](#mandatory-content) and [User-defined columns](#user-defined-columns).

The command requires that the project's dimensions index is up-to-date with the VCF files in the project, and that is why is sort under `divbase-cli dimensions` in the CLI command tree. If you are unsure if the dimensions index is up-to-date, just run `divbase-cli dimensions update` and wait until that job has completed by checking `divbase-cli task-history user`.

The validation command will fetch all sample names from the project dimensions index from the DivBase server and use that to validate that the sample names in the TSV are correct. Misspelled, missing, or otherwise incorrect sample names in the TSV will result in erroneus or even misleading query results, and the validator will help with spotting that.

The following will return errors. These must be fixed if the sidecar TSV should be used in DivBase queries:

- Header formatting: Header row is missing or first column is not #Sample_ID, Duplicate or empty column names

- Tab separation: Row has the wrong number of columns

- `Sample_ID` : Empty Sample_ID,Sample_ID contains a semicolon,Duplicate Sample_ID

- Unsupported characters: no commas in cell values; no hyphens in numerical columns

- Type consistency (numeric and string values): no Mixed types in a column  or in a cell  in a cell (e.g., 1;abc)

- All samples listed in in TSV must exist in the dimensions index

The validator will also raise Warnings. DivBase queries can still be run with these, but the user should review them, and possible address them if so desired:

- Cell value has leading or trailing whitespace (will be stripped by server)

- Samples in the project’s dimensions index not found in the TSV. These samples will not be considered in queries, and that might in fact be what the user wants, espcially if using multiple TSVs. Just be sure to be careful when using this since it will affect the results.

## Query Syntax for sidecar metadata

- TODO: explain warnings, these should be the same as the validator, but this needs to be checked
- TODO: explain when empty results or all results are returned

### Overview: querys are applied as filters on columns in the TSV

Queries on the sidecar sample metadata TSV can be done with the `divbase-cli query tsv` command. The filters that the user want to query on needs entered as a string (i.e. enclosed in quotes, `""`).

The TSV query syntax is `"Key1:Value1,Value2;Key2:Value3,Value4"`, where `Key1:`...`Key2:` are the column header names in the TSV, and `Value1`...`Value4` are the values. Multiple filter values for a key are separated by commas, and multiple keys are separated by semicolons. There can be any number keys and values to filter on, but it is up to the user to write queries that return useful results.
It is possible to exclude a value by prefixing it with a `!` (NOT) operator: `"Key:!Value"`. When mixing inclusive and exclusive filters (e.g. `"Key1:Value1,Value2; Key2:!Value3"`), only the rows that match the positive filters and do not match any of the excluded values will be returned. This can be used to write complex queries.

!!! note
    Please note that semicolon (`;`) is used for different purposes in the TSV (multi-value cells) and in the query syntax (perform queries on multiple columns)!

    Also note that commas are allowed in the query syntax, but are not allowed in the cells in the TSV.

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

The `!`(NOT) operator can be used to exclude specific cell values from a column. When a `!` is used on its own, such as in the command:

```bash
divbase-cli query tsv "Area:!North"
```

it will return all rows that do not contain `North` in the `Area`. Multi-column values that contain `North`, such as a row with e.g. `North;South` will also be excluded by this query.

Note that when inclusive and exclusive are combined (e.g. `"Area:East,!South"`), only rows that match both filters (include `East`, exclude `South`) will be returned in the results.

### Filtering on numerical columns

For numerical columns, it is possible to filter on the following operations:

- **Inequalities**
  Examples: `"Weight:>25"` or `"Weight:>=20,<=40"` or `"Weight:<100"`
  Note: The inequality operator must be expressed relative to the key, i.e. `"Weight:>25"`. The reverse notation `"Weight:25<"` is not supported.
  The syntax only accepts `<=` and `>=` since this is the syntax of Python. The forms `=<` and `=>` are not accepted and will return an error.

- **Range (inclusive)**
  Example: `"Weight:20-40"`

- **Discrete values**
  Example: `"Weight:25,30,35"`

Furthermore, it is possible to combine filters on inequalities, ranges, and discrete values using inclusive OR logic. This means that if any one of the specified conditions is satisfied for a cell, the row will be included in the results. For example:

- `"Weight:<2,4"` returns rows where the value is less than 2 **or** equal to 4
- `"Weight:1-2,4"` returns rows where the value is in the range 1–2 **or** equal to 4
- `"Weight:>5,1-2,4"` returns rows where the value is greater than 5 **or** in the range 1–2 **or** equal to 4
- `"Weight:>10,<2,5-7"` returns rows where the value is greater than 10 **or** less than 2 **or** in the range 5–7

The `!` (NOT) operator can really come to good use for numerical filters:

- `"Weight:!25"` returns rows where the value is not 25.
- `"Weight:>5,!10-15"`  returns rows where the value is greater than 5, but not in the range 10–15.
- `"Weight:!1-2,4"`  returns rows where the value is not in the range 1–2, or is 4.

## Examples of complex queries

Assuming that the sidecar metadata TSV file looks like in the [Example](#example) above, a query like will:

```bash
divbase-cli query tsv "Area:North,West,!South;Weight:>10,<=20,!15,18-22"
```

- include rows where the `Area` column contains either `North` or `West` (also applied to semicolon-separated multi-value cells), **but excludes** any row where `South` is present in the `Area` column—even if `North` or `West` is also present.

- include rows where the `Weight` column is greater than 10, **or** less than or equal to 20, **or** in the range 18–22 (inclusive), **but excludes** any row where Weight is exactly 15 **or** any value in the range 18–22.

There are three samples (rows) that fulfill this, and this is what the query results will return: `S1`, `S4`, and `S5`.

TODOs:

- [TO BE IMPLEMENTED] consider changing the mandatory column name from `Sample_ID` to `Sample`
- [TO BE IMPLEMENTED] what happens if a TSV does not contain all the samples in the DivBase project? There should probably be a warning, but not an error?
- [TO BE IMPLEMENTED] what happens if a sample name is misspelled in the TSV? a warning? can this be checked against the dimensions show?
- [TO BE IMPLEMENTED] what happens if a sample is duplicated in the file. what happens if the sample name is duplicated but not the values (diverging duplicate)?

- [TO BE IMPLEMENTED] what to do if a query references a column that does not exist. E.g. `divbase-cli query tsv "Area:Northern Portugal"` when Area does not exist? This should probably give a warning and not just return nothing

- [TO BE IMPLEMENTED] what to do if a query references a column value. E.g. `divbase-cli query tsv "Area:Northern Portugal"` when Northern Portugal does not exist in the column? This should probably also give a warning and not just return nothing, but nothing is a result here and not a syntax problem...
