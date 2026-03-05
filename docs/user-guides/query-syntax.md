# DivBase Query Syntax for VCF data

TODO

TODO this note is from the dimensions page, but it good to have here too:

!!! Note
    When you submit a query, DivBase will use the state of the VCF Dimensions and the VCF files at that very point in time to produce the query results. It is therefore fine if you or another project member uploads new VCF files to the project while a query is queued or running.

## combined sample metadata and VCF queries

TODO - there is a link to here from the sample metadata guide, so the combined queries should be described in detail here

It is also possible to run a sidecar sample metadata query as part of a VCF query by adding the query as a sting to the flag `--tsv-filter`:

```bash
divbase-cli query bcftools-pipe --tsv-filter "Area:North" --command "view -s SAMPLES"
```
