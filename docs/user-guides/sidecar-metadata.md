# Sidecar Metadata TSV files: creating and querying sample metadata files

DivBase supports that users supply a sidecar TSV (tab separated variables) file with metadata on the samples contained within the VCF files in the DivBase project.

While there are ways for sample metadata to be stored in the VCF itself (see [The Variant Call Format Specification](https://samtools.github.io/hts-specs/VCFv4.5.pdf)), it is not really standardized. Metadata can for instance be specified instance in a global `##SAMPLE` header (once per sample) or in a custom per-variant genotype `FORMAT` field in each variant and sample. The downside of the former is that common tools like `bcftools view` do not filter on the headers; the downside of the latter is that writing the metadata once per variant will result in a lot of repeated data, which in turn leads to elevated file size and processing times as the VCF file scales.

DivBase takes a different approach by decoupling the sample metadata from the VCF data by storing it in a sidecar file. The sidecar TSV can be queried on its own, or together with the VCF files in the DivBase project. The TSV is lightweight and highly extendable (essentially a plain-text form of a spreadsheet). This approach avoids having to read, write, and rewrite metadata to the VCF files and therefore keeps the resource overhead low for the sample metadata.

To be able to accomodate metadata needs for any research project that deals with VCF files, the sidecar sample metadata TSV and filtering in DivBase has been designed to be very open-ended and user-defined. As long as a few format and filter syntax requirements are followed, the user is free to design their metadata TSV as they like using the format: column names represent metadata categories and rows represent the samples found in the VCF files in the DivBase project. However, this flexibility put the responsibility on the user that spelling and values in columns and rows are correct: if not, the sample metadata filters will return incomplete or unintended results.

!!! Notes
    There is a CLI command to help check that a user-defined sample metadata TSV file aligns with the requirements described on this page. This validator tool will be [described in its own section below](#validating-a-sidecar-metadata-tsv-with-divbase-cli), but, in short, it can be run with:

    ```bash
    divbase-cli dimensions validate-metadata-file path/to/your/sample_metadata.tsv
    ```

This guide contains sections on how to [Create a sample metadata TSV](#creating-a-sidecar-sample-metadata-tsv-for-a-divbase-project), and [How to run queries on sample metadata TSV files](#query-syntax-for-sidecar-metadata). Instructions on how to run combined sample metadata and VCF data queries are found on the separate page on [DivBase Query Syntax for VCF data](query-syntax.md).

!!! Warning
    All instructions regarding running DivBase queries, generating sample metadata templates, and validating sample metadata TSV files required that the project's VCF dimensions index is updated against the current versions of the VCF files in the project's data store. This can be assured by running the command:

    ```bash
    divbase-cli dimensions update
    ```

    Depending on the number and sizes of the VCF files, this can take a little time. To check the status of the dimensions update job, use the command:

    ```bash
    divbase-cli task-history user
    ```

## Creating a sidecar sample metadata TSV for a DivBase project

If the dimensions VCF files in the project have been cached in DivBase, a template metadata file with the sample names pre-filled can be created with:

```bash
divbase-cli dimensions create-metadata-template
```

!!! Note
    There can be multiple TSVs in the same project and it is possible to call them for the queries with the `--metadata-tsv-name` flag. If not specified, it the default `sample_metadata.tsv` will be assumed. It is up to the user if the want to have multiple TSVs in the same project to organise their metadata in a specific way. It is allowed to have duplicate sample names and metadata across multiple TSV files, since only one TSV can be called per query. It is recommended to a have a master TSV that contains all samples from all the VCFs in the project: querying on TSVs that contain subsets of all sample names is possible, but will sample names not included in the TSV used for the query will be disregarded for the query.

### Sidecar TSV format requirements

To be able to accommodate a variety of metadata needs, DivBase does not enforce a strict schema for the sidecar sample metadata TSV file since it is designed to contain user-defined columns. Instead, there are a few mandatory requirements and some best-practices for defining columns.

#### Mandatory content

1. The first row must be a header row, start with `#`, and the first column must be named `Sample_ID`.
2. The `Sample_ID` column must contain the exact names of the samples as they are spelled in the VCF files. Sample names need to occur uniqely in the TSV: only one row per sample name in the `Sample_ID` column, no duplicates allowed. This will already be handled if user has run a `divbase-cli dimensions update` job and, after its completion, has generated a pre-filled template with: `divbase-cli dimensions create-metadata-template`
3. The `Sample_ID` column can only contain one sample name per row. This is different from the user-defined columns that can take arrays of values for each cell in a column using semicolons (`;`) as delimters. `Sample_ID` values can also not be empty.

4. Every column need to be tab separated for all rows, including the header.

#### User-defined columns

After the `Sample_ID` column has been populated, users can add any columns and values to the TSV.

!!! Warning
    It is the user's responsibility to ensure that the spelling of column headers and values is consistent. When filtering on the sidecar metadata, the exact spelling of a column name must be used for the filters. This includes matching upper and lower case letters.

To ensure that user-defined metadata can be used in DivBase, we ask you follow the following constraints and considerations:

1. The user-defined columns can be **either** numeric **or** string type. A column is classified as numeric only if all values can be parsed as numbers (including individual parts in semicolon-separated cells). If any value in a column is non-numeric, the entire column is treated as a string column. This means a column with values like "8", "1a", "5a" will be treated as string column even though some values look numeric. The DivBase backend uses [`Pandas`](https://pandas.pydata.org/) to automatically infer column type based on its data, so there is no need to specify in the TSV whether the values are numerical or string.
2. Semicolon-separated values are supported in TSV cells to represent arrays of values. This allows users to have samples that can belong to multiple values in the same column. For instance belong to two different groups or categories. This works with both numerical and string data (e.g. "2;4;21" or "North; North-West"). Note that this might make the process of writing queries more complex than if just a single value is used for each cell. **Important:** Semicolons (`;`) are the only supported delimiter for multi-value cells. DivBase uses commas (`,`) in the [Query syntax](#query-syntax-for-sidecar-metadata) for a different purpose (separating filter values in queries).
3. Special characters like hyphens (`-`) and commas (`,`) are allowed, but will cause the column to be treated as a string column. String columns cannot be filtered using numeric operators (see details in [Filtering on numerical columns](#filtering-on-numerical-columns)) and will raise warnings. For example, values like "1-2" or "1,2" will be interpreted as strings, not numeric ranges or multi-value fields. If you intend to store multiple numeric values in a cell, use semicolons (e.g., "1;2"). For decimals, use English decimal notation with a period (e.g., "3.14") and not a comma.
4. The only characters with special structural meaning in DivBase sidecar metadata TSV files are `#` (for header comments), `;` (for multi-value cell separation), and `\t` (tab, for column separation). Other special characters are generally supported in data values, but be aware that Your Mileage May Vary. Some common cases that have been tested and are supported include diacritic unicode characters like `å`, `ä`, `ö`, and hyphens in strings (e.g., `North-West`).
5. Leading and trailing whitespaces are removed by the DivBase backend in order to ensure robust filtering and pattern matching. Whitespaces inside strings will be preserved. For instance: " Sample 1 " will be processed as "Sample 1".

#### Example

This example illustrates how a sidecar sample metadata TSV can look like. The mandatory requirement are fulfilled (heading with `#`, `Sample_ID` column, tab-separated file). The user-defined column contains examples of a numerical column (`Population`) and a string column (`Area`). In some cells, semicolons (`;`) are used to assign multiple values to the same sample and column.

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

The validation runs on the user's local computer and not as a job on the DivBase server. It is intended to be used on sidecar metadata TSV files before they are uploaded to the DivBase project. The validator will check the formatting requirements as described in [Mandatory contents](#mandatory-content) and [User-defined columns](#user-defined-columns).

The command requires that the project's dimensions index is up-to-date with the VCF files in the project, and that is why it is sorted under `divbase-cli dimensions` in the CLI command tree. If you are unsure if the dimensions index is up-to-date, just run `divbase-cli dimensions update` and wait until that job has completed by checking `divbase-cli task-history user`.

The validation command will fetch all sample names from the project dimensions index from the DivBase server and use that to validate that the sample names in the TSV are correct. Misspelled, missing, or otherwise incorrect sample names in the TSV will result in erroneous or even misleading query results, and the validator will help with spotting that. Several of the checks that the validator performs are also done at the start of a sample metadata query, but this sample name check is currently only done by the validator.

The following will return **Errors**. These must be fixed for the sidecar TSV be used with DivBase queries:

- File not found, unreadable, or empty: If the TSV file path is missing, misspelled, or the file cannot be opened, validation will fail. Empty TSV files are also not allowed.

- Header formatting: Header row is missing or first column is not `#Sample_ID`, duplicate or empty column names.

- Tab separation: Row has the wrong number of columns. (Note: This check is only done in the validator! It is currently not part of the checks at the start of a sample metadata query.)

- `Sample_ID` column issues: Empty value, value containing a semicolon, rows with duplicate sample names.

- Samples in TSV not found in project dimensions index: All samples listed in the TSV must exist in the project's dimensions index. If a sample is known to be in a VCF file in the DivBase project but is missing from the VCF dimensions index, the user needs to run `divbase-cli dimensions update` to submit an update job and then try the validator again after the job has finished.

!!! Note
    The formatting errors listed above are also enforced by the DivBase query engine when loading the metadata file for queries (except checking tab separation which is a validator-specific check). This means that even if the validator is not run before upload, the query engine will analyse the file content and report issues as errors. Detected Errors are different from Warnings in that errors will result in queries not even being run.

The validator will also raise **Warnings**. DivBase queries can still be run with Warnings, but the user should review them, and possible address them if so desired:

- Cell value has leading or trailing whitespace (will be stripped by DivBase when a query is run)

- Samples in the project's dimensions index not found in the TSV. These samples will not be considered in queries, and that might in fact be what the user wants, especially if using multiple TSVs. Just be sure to be careful when using this since it will affect the results.

- Mixed-type columns (a column with numeric and string values, e.g. "8", "1a", "5a") and semicolon-separated cells with mixed types (e.g., "1;abc"). They are allowed but the user should keep in mind that since they will be treated as string columns, numeric query operations (ranges, inequalities) will not work on these columns.

- Hyphens in values that look like range notation (e.g., "1-2") in columns that otherwise contain numeric values. The same goes for commas (e.g. "1,2"). The warning message will ask the user if they intended this to be a multicolumn value which should use semicolons as delimiters.

## Query Syntax for sidecar metadata

This section describes how to query on the sample metadata file itself. The same syntax used here will also be used when running combined sample metadata and VCF data queries; how to do that is covered in [DivBase Query Syntax for VCF data](query-syntax.md).

### Overview: querys are applied as filters on columns in the TSV

Queries on the sidecar sample metadata TSV can be done with the `divbase-cli query tsv` command. The filters that the user wants to query on need to be entered as a string (i.e. enclosed in quotes, `""`).

The TSV query syntax is `"Key1:Value1,Value2;Key2:Value3,Value4"`, where `Key1:`...`Key2:` are the column header names in the TSV, and `Value1`...`Value4` are the values. Multiple filter values for a key are separated by commas, and multiple keys are separated by semicolons. There can be any number keys and values to filter on, but it is up to the user to write queries that return useful results.
It is possible to exclude a value by prefixing it with a `!` (NOT) operator: `"Key:!Value"`. When mixing inclusive and exclusive filters (e.g. `"Key1:Value1,Value2; Key2:!Value3"`), only the rows that match the positive filters and do not match any of the excluded values will be returned. This can be used to write complex queries.

!!! note
    Please note that semicolon (`;`) is used for different purposes in the TSV (for denoting multi-value cells) and in the query syntax (for performing queries on multiple columns)!

For example, if the user wants to query the TSV on column `Area` for all samples that contain the value `North`,:

```bash
divbase-cli query tsv "Area:North"
```

Please also see the documentation on [DivBase Query Syntax for VCF data](query-syntax.md) for more details on how that command works.

Filtering is inclusive by default. This applies both for the filter values and the cell values:

- If a filter contains multiple values, e.g. `"Area:North,West"`, the row is included if at least one of the filter values matches any value in the cell. I.e. a row with `North`, and a row with `West` will both be returned from this filter.
- If a cell in the TSV contains multiple values separated by a semicolon as explained in [User-defined columns](#user-defined-columns) (e.g., `North;West`), the row is included if any of those values match the filter. Filters with `"Area:North"`, `"Area:West"`, and `"Area:North,West"` will all return the row with the array value `North;West`.

!!! note
    To reiterate what was written in the [User-defined columns](#user-defined-columns) section above: it is the user's responsibility to ensure that the spelling of column headers and values is consistent. When filtering on the sidecar metadata, the exact spelling must be used for the filters.

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

A TSV column is considered as numeric in DivBase only if all cell values — including each individual part within semicolon-separated cells (e.g. `1;3;5`) — can be parsed as a number. For example:

- A column with values `1`, `2;4`, `3`, `1;3;5` is considered numeric since all elements are numbers. All numeric operations below (inequalities, ranges, discrete) are fully supported on this column.

- A column with values `1;1-2`, `3`, `5` is considered a string column since the part `1-2` cannot be parsed as a number. Only exact string matching is supported for this column.

- A column with values `8`, `1a`, `5a`is considered a string column since it has mixed types (`8` is numeric, the others are strings). Only exact string matching is supported for this column.

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

!!! Tip
    Numeric operations such as inequalities like `>25`, and ranges like `20-40` are fully supported for semicolon-separated numeric columns as long as every semicolon separated part (`part;part`) in every cell in the column is a valid number. For instance: a `Population` column with values `1`, `2;4`, `1;3;5`; in this case a query like `divbase-cli query tsv "Population:>3"` will correctly match cells like `2;4` and `1;3;5`.

### Query Warnings: spotting potential issues with the TSV or the query filter

When running a sample metadata query in DivBase, the system will check the TSV and the query filter for the constraints and considerations described throughout this guide. If errors are encountered, the query will not run and a message with details on what went wrong will be return to the user. Warnings, however, will not stop queries from running, but indicate that the user should carefully review the results.

Reviewing the Warnings to judge if they are relevant or not is essential to avoid unintended query results. The following are treated as Warnings by DivBase queries (and by the TSV validator).

!!! note
    Note that the warnings will only occur when query contains a combination of filter and column values that warrants a warning: filtering only on Column A will not check for warnings in Column B. Use the [TSV validator](#validating-a-sidecar-metadata-tsv-with-divbase-cli) if you want to check for issues across all columns.

- **Comparison operators on string/mixed-type columns**: DivBase comparison operators (`>`, `<`, `>=`, `<=`) only work on numeric columns. If you use them on a string or mixed-type column — whether with a numeric operand (e.g., `Population:>5`) or a string operand (e.g., `Area:>North`) — DivBase will warn that comparison operators are not supported on string columns. Use exact string matching instead (e.g., `Area:North` or `Population:8,1a`).

- **Mixed-type column information**:
When filtering on a mixed-type column with valid string matching, DivBase will inform you that the column is treated as string and comparison operators are not available. This is mainly to make the user aware of this.

- **Column not found**:
If the filter references a column that does not exist in the TSV, DivBase will warn and skip that filter condition.

- **No matching values**:
If none of the filter values match any values in the column, DivBase print a warning. This can indicate a typo in the filter value, or just that the specific filter combination filtered away all samples..

### Examples of complex queries

Assuming that the sidecar metadata TSV file looks like in the [Example](#example) above, a query like will:

```bash
divbase-cli query tsv "Area:North,West,!South;Weight:>10,<=20,!15,18-22"
```

- include rows where the `Area` column contains either `North` or `West` (also applied to semicolon-separated multi-value cells), **but excludes** any row where `South` is present in the `Area` column—even if `North` or `West` is also present.

- include rows where the `Weight` column is greater than 10, **or** less than or equal to 20, **or** in the range 18–22 (inclusive), **but excludes** any row where Weight is exactly 15 **or** any value in the range 18–22.

There are three samples (rows) that fulfill this, and this is what the query results will return: `S1`, `S4`, and `S5`.
