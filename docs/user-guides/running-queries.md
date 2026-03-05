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

To ensure that users get good performance out of DivBase, technical metadata from all VCF files that are uploaded to a DivBase project will be indexed on the DivBase server. For each VCF file name, the server will extract the number and name of all samples, the number and name of all scaffolds, and the number of variants. This technical VCF metadata is referred to as VCF Dimensions in DivBase. The full details of how to use this are described on the page: [VCF Dimensions](vcf-dimensions.md).

DivBase will use the VCF Dimensions whenever a user submits a query to the project. Instead of transferring VCF files internally and reading them each time a user sends a query to DivBase, the server will instead use the VCF Dimensions cache, which is a much quicker process.

The VCF Dimensions must be updated each time a new VCF file or a new version of a VCF file is uploaded to the project. This needs to be done manually by running the command:

```bash
divbase-cli dimensions update
```

This command is sent to the DivBase job queue system, since this process can potentially take some time for larger VCF files. DivBase uses `bcftools` to extract the VCF dimensions metadata, which is a fast software. Nevertheless, the time it will take to extract the VCF dimensions will scale proportionally to the VCF file size. You can see this as a once-per-VCF-file-version investment: it will take a little time here up-front, but all subsequent queries to this file will be much faster than they otherwise would be.

Please wait until the VCF dimensions job has finished for your project before you make any queries to DivBase. You can check the status of the job with:

```bash
divbase-cli task-history user
```

To see the current VCF dimensions for a project, use:

```bash
divbase-cli dimensions show
```

This will print the currently cached VCF Dimensions in DivBase per VCF file name ([see example here](vcf-dimensions.md#dimensions-show))

The VCF Dimensions of the project should now be up-to-date with the VCF files in the projects data store. You are now ready to run queries.

!!! Note
    Remember that if you or another project member uploads a VCF file to the project, you will need to run `divbase-cli dimensions update` again before running queries.

    When you submit a query, DivBase will use the state of the VCF Dimensions and the VCF files at that very point in time to produce the query results. It is therefore fine if you or another project member uploads new VCF files to the project while a query is queued or running.

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
