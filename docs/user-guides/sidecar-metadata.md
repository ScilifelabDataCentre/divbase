# Sidecar Metadata TSV files: creating and querying sample metadata files

TODO

- rationale: what is this and how can it be used

## Creating a sidecar TSV for a DivBase project

If the dimensions VCF files in the project have been cached in DivBase, a template metadata file with the sample names pre-filled can be created with:

```bash
divbase-cli dimensions create-metadata-template
```

Note! there can be multiple TSVs in the same project and it is possible to call them for the queries with the `--metadata-tsv-name` flag.

TODOs:

- [TO BE IMPLEMENTED] consider changing the mandatory column name from `Sample_ID` to `Sample`
- [TO BE IMPLEMENTED] what happens if a TSV does not contain all the samples in the DivBase project? There should probably be a warning, but not an error?
- [TO BE IMPLEMENTED] what happens if a sample name is misspelled in the TSV? a warning? can this be checked against the dimensions show?

- [TO BE IMPLEMENTED] what happens if a sample is duplicated in the file. what happens if the sample name is duplicated but not the values (diverging duplicate)?

1 mandatory column: sample name.

any other columns are optional and user defined

it is possible to have more than one sidecar sample metadata file in each DivBase project.

example

```
TODO add example here
```

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

TODO write pytests that ensure that these numerical filters work

## Trying out a query

```bash
divbase-cli query tsv "Area:Northern Portugal"
```

- [TO BE IMPLEMENTED] what to do if a query references a column that does not exist. E.g. `divbase-cli query tsv "Area:Northern Portugal"` when Area does not exist? This should probably give a warning and not just return nothing

- [TO BE IMPLEMENTED] what to do if a query references a column value. E.g. `divbase-cli query tsv "Area:Northern Portugal"` when Northern Portugal does not exist in the column? This should probably also give a warning and not just return nothing, but nothing is a result here and not a syntax problem...
