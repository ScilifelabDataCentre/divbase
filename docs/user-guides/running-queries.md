# Running Queries: Overview

DivBase allows users to submit queries to checkout sample metadata and VCF data from the DivBase project data store. There types of queries in DivBase are:

- Sample metadata query
- VCF data query
- Combined sample metadata and VCF data query

The sample metadata is a user defined, sidecar TSV (tab separated variables) file where users can add any custom metadata as columns. The VCF queries use `bcftools` to subset data from the files in the data store based on user defined filters. The VCF queries can be combined with a sample metadata query to allow users to filter both on data in the sidecar metadata and VCF files in a single query.

TODO

- [TO BE IMPLEMENTED] think about the `divbase-cli query` commands. it would make sense to use `sample-metadata` (perhaps with `tsv` as a short form) and `vcf`

## Side car sample metadata queries

For more details, see [Sidecar Metadata TSV files: creating and querying sample metadata files](sidecar-metadata.md).

how to write the sample metadata TSV (+template)

how to use dimensions show to get all samples in the project

how to write query

## Before any VCF queries: run the VCF dimensions command

For more details see [VCF Dimensions caching](vcf-dimensions.md)

For performance reasons and to ensure query feasibility, key metadata from the VCF files must first be cached in DivBase.
For each VCF in the DivBase project, the system extracts file name, number and name of all samples, number and name of all scaffolds, number of variants,

DivBase will use this whenever a user submits a query to the project. For instance, the user might make a sample metadata filtering that results in only certain samples. The system knows which file names each requested sample are located in, and will ensure that only those files will be transferred to the worker.

example
dimensions show

## VCF queries

See also [DivBase Query Syntax for VCF data](query-syntax.md), [How to create efficient DivBase queries](how-to-create-efficient-divbase-queries.md), and [Tutorial: Running a query on a public dataset](tutorial-query-on-public-data.md).

can be run with or without sample metadata filtering

for sample metadata linked VCF queries, it can be good to do a dry run first [TO BE IMPLEMENTED, at the moment it needs to be run seperatelly]

how to write a query

how to optimize a query (see separate markdown)
