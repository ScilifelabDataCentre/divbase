# Running Queries: Overview

DivBase allows users to submit queries to checkout sample metadata and VCF data from the DivBase project data store. The three types of queries in DivBase are:

- Sample metadata query
- VCF data query
- Combined sample metadata and VCF data query

The sample metadata is a user defined, sidecar TSV (tab separated variables) file where users can add any custom metadata as columns. The VCF queries use `bcftools` to subset data from the files in the data store based on user defined filters. The VCF queries can be combined with a sample metadata query to allow users to filter both on data in the sidecar metadata and VCF files in a single query.

A few quick things that are good to know when it comes to running DivBase queries:

- The queries are run as jobs in the DivBase job queue management system since they can potentially take some time to run, depending on the VCF file size and the query itself.

- The system will use the latest version the files for all queries.

- The VCF Dimensions for the project needs to be up-to-date before running any queries, as described in [Before running any queries: run the VCF dimensions command](#before-running-any-queries-run-the-vcf-dimensions-command) below.

This page is meant as medium detailed description of these query commands, and therefore sits in-between the [Quick Start Guide](quick-start.md) and the full guides found here in the Running Queries chapter in terms of level-of-detail. The sections below contain links to the full guides.

## Before running any queries: run the VCF dimensions command

TODO - finish writing this section

For more details see [VCF Dimensions caching](vcf-dimensions.md)

For performance reasons and to ensure query feasibility, key metadata from the VCF files must first be cached in DivBase.
For each VCF in the DivBase project, the system extracts file name, number and name of all samples, number and name of all scaffolds, number of variants,

DivBase will use this whenever a user submits a query to the project. For instance, the user might make a sample metadata filtering that results in only certain samples. The system knows which file names each requested sample are located in, and will ensure that only those files will be transferred to the worker.

example
dimensions show

how to use dimensions show to get all samples in the project

## Sidecar sample metadata queries

DivBase supports that users store extensive sample metadata in a separate TSV file. This metadata can be queried on its own or in combination with VCF data. Users are free to define their own metadata as they see if: column names represent metadata categories and rows represent the samples found in the VCF files in the DivBase project. The TSVs need to follow a few mandatory requirements, but no strict metadata schema is enforced, This allows DivBase to accomodate a variety of research projects with different metadata needs.

The metadata can be queried on its own to learn which samples that fulfil a certain metadata query and the VCF files the samples are present in. The same query syntax is used in the combined sample metadata and VCF data queries, and user can use the dedicated sample metadata query command as a dry-run before running full combined query to ensure that the metadata query produces the results the user intended.

For instructions on how to create the sidecar sample metadata TSV files and how to run sample metadata queries, see the guide on [Sidecar Metadata TSV files: creating and querying sample metadata files](sidecar-metadata.md). The guide also describes the CLI commands that specifically relate to the sample metadata TSV files. These are:

```bash
divbase-cli dimensions create-metadata-template

divbase-cli dimensions validate-metadata-file path/to/your/sample_metadata.tsv

divbase-cli files upload path/to/your/sample_metadata.tsv

divbase-cli query tsv "Area:Northern Portugal"
```

## VCF queries

TODO - finish writing this section

Run the query on all VCF files in the DivBase project unless specified. There are two ways to specify the files: either as direct input to the `bcftools` command, or by combining the VCF query with a sample metadata query to determine which VCF files to use.

See also [DivBase Query Syntax for VCF data](query-syntax.md), [How to create efficient DivBase queries](how-to-create-efficient-divbase-queries.md), and [Tutorial: Running a query on a public dataset](tutorial-query-on-public-data.md).

can be run with or without sample metadata filtering.

for sample metadata linked VCF queries, it can be good to do a dry run first [TO BE IMPLEMENTED, at the moment it needs to be run seperatelly]

how to write a query

how to optimize a query (see separate markdown)

TODO COPIED OVER FROM QUICKSTART REIMPLEMENT

DivBase only allows `bcftools view` in its query syntax and no other `bcftools` commands. The `merge`, `concat`, and `annotate` commands are used when processing a query, but should not be defined by the user.

## Combined sample metadata and VCF data query

TODO - finish writing this section

Uses a sample metadata query to identify the VCF files in the DivBase project to run the VCF queries on.
