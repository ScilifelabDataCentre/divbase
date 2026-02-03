# Sidecar Metadata TSV files: creating and querying sample metadata files

TODO

- rationale: what is this and how can it be used

## Creating a sidecar TSV for a DivBase project

Note! there can be multiple TSVs in the same project and it is possible to call them for the queries with the `--metadata-tsv-name` flag.

TODOs:

- [TO BE IMPLEMENTED] consider changing the mandatory column name from `Sample_ID` to `Sample`
- [TO BE IMPLEMENTED] CLI command to generate template (empty template and template with the samples from the DivBase project pre-filled). Pre-filling the template will require that dimensions update has been run
- [TO BE IMPLEMENTED] what happens if a TSV does not contain all the samples in the DivBase project? There should probably be a warning, but not an error?
- [TO BE IMPLEMENTED] what happens if a sample name is misspelled in the TSV? a warning? can this be checked against the dimensions show?

1 mandatory column: sample name.

any other columns are optional and user defined

example

```
TODO add example here
```

## Query Syntax for sidecar metadata

from `divbase-cli query tsv -h` docstring:
String consisting of keys:values in the tsv file to filter on. The syntax is 'Key1:Value1,Value2;Key2:Value3,Value4', where the key are the column header names in the tsv, and values are the column values. Multiple values for a key are separated by commas, and multiple keys are separated by semicolons. When multple keys are provided, an intersect query will be performed. E.g. 'Area:West of Ireland,Northern Portugal;Sex:F'.

- [TO BE IMPLEMENTED] filtering based on ranges (ints and maybe floats) and not just on strings, e.g. in range 31 - 50. etc...
- [TO BE IMPLEMENTED] add more Set Operations: union, intersection, difference, symmetric difference. Be clear on the default behaviour

## Trying out a query

```bash
divbase-cli query tsv "Area:Northern Portugal"
```

- [TO BE IMPLEMENTED] what to do if a query references a column that does not exist. E.g. `divbase-cli query tsv "Area:Northern Portugal"` when Area does not exist? This should probably give a warning and not just return nothing

- [TO BE IMPLEMENTED] what to do if a query references a column value. E.g. `divbase-cli query tsv "Area:Northern Portugal"` when Northern Portugal does not exist in the column? This should probably also give a warning and not just return nothing, but nothing is a result here and not a syntax problem...
