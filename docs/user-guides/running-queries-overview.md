# Running Queries: Overview

A core functionality of DivBase is to allow users to checkout user-defined subsets of the VCF data and sample metadata contained in their project. In DivBase, the action of checking out data like this is called a _query_.

There are three query workflows:

- Query sidecar metadata only
- Query VCF data only
- Combine sidecar metadata filtering with a VCF query

All queries depend on that the prerequisites are fulfilled, as described in the next section.

## Prerequisites

First, ensure that the VCF files fulfill the requirements described in [Working with VCF Files in DivBase](vcf-files.md) and have been uploaded to the DivBase project. Then ensure that the project's VCF dimensions cache is up-to-date by running:

```bash
# Submit a job to that ensure that the projecs' dimensions cache is up-to-date
divbase-cli dimensions update --project <PROJECT_NAME>

# To check the progress of the dimensions update job, use
divbase-cli task-history user
```

To be able to query on sample metadata, you also need to upload a TSV file to the project, as described in the guide on [Sidecar Metadata TSV files](sidecar-metadata.md).

Once the VCF dimensions update job is finished and the dimensions cache of the project is up-to-date, continue with one of the query paths below.

## Query paths

| What you want to do | Command | Full guide |
|---|---|---|
| Find samples/files from metadata TSV | `divbase-cli query tsv "<FILTER>"` | [Sidecar Metadata TSV files](sidecar-metadata.md) |
| Submit a VCF subset job | `divbase-cli query vcf --command "view ..."` | [DivBase VCF query syntax](vcf-query-syntax.md) |
| Use metadata filter + VCF filter together | `divbase-cli query vcf --tsv-filter "<FILTER>" --command "view ..."` | [DivBase VCF query syntax](vcf-query-syntax.md) |

## Minimal examples

```bash
# Metadata query only
divbase-cli query tsv "Area:North"

# VCF query only
divbase-cli query vcf --command "view -r 21:15000000-25000000"

# Combined metadata + VCF query
divbase-cli query vcf \
  --tsv-filter "Area:North,West" \
  --command "view -r 21:15000000-25000000"
```

## What happens after query submission?

- `query tsv` directly returns the sample IDs and VCF filenames in the project matching the query.
- `query vcf` submits an asynchronous job and returns a task ID.
- On successful VCF jobs, DivBase uploads a results VCF file with the subset data in the project's data storage.

Check VCF job status with:

```bash
# By the job ID
divbase-cli task-history id <JOB_ID>

# Or by all the jobs submitted by the user
divbase-cli task-history user
```

## Read next

- [VCF Dimensions caching](vcf-dimensions.md)
- [Sidecar Metadata TSV files: creating and querying sample metadata files](sidecar-metadata.md)
- [DivBase VCF query syntax](vcf-query-syntax.md)
- [How to create efficient DivBase queries](how-to-create-efficient-divbase-queries.md)
- [Tutorial: Running a query on a public dataset](user-guides/tutorial-query-on-public-data.md)
