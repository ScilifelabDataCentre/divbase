# Sidecar sample metadata TSV queries

The query engine for filtering on user-defined sidecar TSVs is centralised around the class `SharedMetadataValidator` in the `divbase-lib` package. This class validated that a sidecar TSV is compatible with DivBase, and is called by the client-side validator and by the server-side query engine (`SidecarQueryManager`). Important to understand is that the `SharedMetadataValidator` was designed to not raise errors or reports warnings: it only collects them. Then the CLI and query engine each acts on the information.

The client side uses the results from `SharedMetadataValidator` to print stats, errors, and warnings in the user's terminal.

The server side uses Enums in the `SharedMetadataValidator` results to raise errors and report warnings to the user. This means that the `SidecarQueryManager` only need to be concerned about acting on the Enum and not the error/warning messages themselves, thus avoiding having to rely on (potentially fragile) string-matching.

The syntax for the TSV and the query fileters is described in the [Sidecar Metadata TSV files: creating and querying sample metadata files](../user-guides/sidecar-metadata.md).

## Script to inspect the dataframe

Despite all the tests, there might be cases where the Pandas dataframe does not look like you'd think after ingesting the TSV. Such issues can be difficult to debug. To help with that this script can be used to calls the `SharedMetadataValidator` dataframe logic for TSV ingestion and print the Pandas dataframe to terminal:

```bash
python scripts/tsv_to_dataframe.py --tsv path/to/metadata.tsv

# Try it with example fixture
python scripts/tsv_to_dataframe.py --tsv tests/fixtures/sample_metadata.tsv
```
