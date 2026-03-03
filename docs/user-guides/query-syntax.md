# DivBase Query Syntax for VCF data

TODO

## combined sample metadata and VCF queries

TODO - there is a link to here from the sample metadata guide, so the combined queries should be described in detail here

It is also possible to run a sidecar sample metadata query as part of a VCF query by adding the query as a sting to the flag `--tsv-filter`:

```bash
divbase-cli query bcftools-pipe --tsv-filter "Area:North" --command "view -s SAMPLES"
```
