# DivBase VCF query syntax

Users can checkout subsets of their VCF data from their DivBase project using the command `divbase-cli query vcf`. This data checkout is run as an asynchronous job that is sent to the queue on the DivBase server, and eventually run once there are idle resources to process the job. Users will get a job ID when they submit their query to the task queue, and can view the status of the job with for instance `divbase-cli task-history id <JOB_ID>`.

The processing of the VCF files on the DivBase server is done with [`bcftools`](https://github.com/samtools/bcftools). DivBase will detect the VCF files in the project's data store that are needed for the query; if more than one VCF file is needed, DivBase will ensure that the files are compatible with each other according to the requirements of `bcftools` and ensure that a single results file with the subset data is returned to the user by running `bcftools merge` and `bcftools concat` on the intermediate files as needed. The result is a single VCF file that is uploaded to the projects data store and named after the job ID.

Users can query the VCF data in their project with or without combining it to a [sample metadata query](docs/user-guides/sidecar-metadata.md).

Example of a VCF query that identifies the samples and VCF files to filter on in the project's datastore and then applies a subset based on genomic range:

```bash
divbase-cli query vcf --tsv-filter  "Area:North,West;Weight:>10" --command "view -r 21:15000000-25000000"

# This will return the job ID of the submitted job. Example:
# Job submitted successfully with task id: 123

# Job status can be viewed with e.g.
divbase-cli task-history id 123

```

The outcome of a DivBase VCF query is a single results file with merged/concatenated data that fulfills all user-defined filters.

!!! Note
    When you submit a query, DivBase will use the state of the VCF Dimensions and the VCF files at that very point in time to produce the query results. It is therefore fine if you or another project member uploads new VCF files to the project while a query is queued or running.

## Prerequisites

Ensure that VCF files and sample metadata TSV file(s) are uploaded to the DivBase project, and that VCF dimensions cache is up-to-date, as described in e.g. the [Running Queries Overview - Prerequisites](running-queries-overview.md#prerequisites).

DivBase uses [`bcftools`](https://github.com/samtools/bcftools) to subset VCF data. The DivBase VCF query syntax is based on `bcftools view`, which is described in the [bcftools manual](https://samtools.github.io/bcftools/bcftools.html#view). If you are not familiar with `bcftools view`, you might want to take some time to study the different options. The commands used for DivBase VCF queries are described in more detail in the [Writing the bcftools command argument](#writing-the-bcftools-command-argument) section below.

TODO bcftools rules, should probably go into vcf-files.md

## `divbase-cli query vcf` command structure

```bash
divbase-cli query vcf \
  [--tsv-filter "<FILTER>"] \
  [--samples "<ID1,ID2,...>"] \
  [--samples-file path/to/samples.txt] \
  --command "<BCFTOOLS_VIEW_PIPE>" \
  [--metadata-tsv-name <FILENAME.tsv>] \
  [--project <PROJECT_NAME>]
```

- Required: `--command`
- Optional: exactly one sample-selection mode (`--tsv-filter` OR `--samples` OR `--samples-file`), or none
- Optional: `--metadata-tsv-name` (mainly needed with `--tsv-filter`)
- Optional: `--project` if user config default exists

Explain mutual exclusivity clearly:

- `--tsv-filter`, `--samples`, and `--samples-file` are mutually exclusive

## Sample and VCF file selection

To run a VCF data query, the user needs to input the samples and/or VCF files to perform the query on, as well as the `bcftools view` command(s) (described in [Writing the bcftools command argument](#writing-the-bcftools-command-argument)) below. DivBase supports different ways to input the samples and/or VCF files:

### Metadata-driven sample and VCF file selection (--tsv-filter)

This is a combined sample metadata and VCF data query, that allows users to let the results of the metadata query (samples and VCF files) be automatically used for a VCF data query. The VCF queries in DivBase was designed with this in mind, since it augments regular `bcftools` subsetting with metadata-guides filtering.

The metadata query is input with `--tsv-filter` argument and uses the TSV format and filter syntax described in the guide on [sample metadata queries](sidecar-metadata.md). They system will by default look for `sample_metadata.tsv` in the project's data storage, but this can be overridden to use another TSV in the data storage using `--metadata-tsv-name <MY_SAMPLE_METADATA.TSV>`.

```bash
divbase-cli query vcf --tsv-filter  "Area:North,West;Weight:>10" --command "view -r 21:15000000-25000000"
```

For example, the system might find that the samples that fulfill the [sample metadata query](docs/user-guides/sidecar-metadata.md) set with `--tsv-filter` are, say, `S2`, `S5`, `S28`, `S108` and that they are described in the files `file1.vcf`, `file3.vcf`, `file4.vcf`. The DivBase server will then act on only these three files and subset based on the four samples.

!!! Tip
    Before using `--tsv-filter` in `query vcf`, you can do a dry-run of metadata query to ensure that the metadata query returns the expected samples and VCF files:

    ```bash
    divbase-cli query tsv "Area:North,West;Weight:>10"
    ```

### Sample selection from direct input (--samples)

Users can also run VCF data query without metadata queries by defining which samples and/of VCF files to subset on. To get an overview of the VCF files and samples in the project, ensure that the [VCF dimensions cache is up-to-date](running-queries-overview.md#prerequisites) and run the following:

```bash
divbase-cli dimensions show
```

To just get all samples that are available for a project, use `divbase-cli dimensions show --unique-samples`

TODO add dimensions show command to get all VCF files

To specify the samples on the command line, use the option `--samples`:

```bash
divbase-cli query vcf --samples "S1,S2,S10,S239" --command "view -r 21:15000000-25000000"
```

As long as the samples are present in the DivBase project, the server will automatically find out the VCF files it needs to process for the VCF query, by reading the VCF dimensions cache.

### Sample selection from file (--samples-file)

An alternative to `--samples` for non-metadata driven VCF queries is to provide a file that contain all the sample names to be used for the query. This can be convienent for queries with a higher number of samples. For example:

```bash
divbase-cli query vcf --samples-file samples_for_my_query.txt --command "view -r 21:15000000-25000000"
```

The file needs to be a plain file, such a `.txt` and contain one sample per row. Example:

```text
# samples_for_my_query.txt
S1
S2
S10
S239
```

TODO ensure that there can be comments in the plain files used for `--samples-file`

### VCF file selection from direct filenames

TODO not implemented

could look like:

```bash
divbase-cli query vcf --vcf-files "file1.vcf, file5.vcf" --command "view -r 21:15000000-25000000"
```

if this is combined with `--samples "S1,S2,S10,S239" --vcf-files "file1.vcf, file5.vcf"` there might be many queries without results. Support for this would need to be considered in the future.

### No selection flags (all samples and files)

TODO: we might not want to support this, it would be easier for the user to download all the files and merge them themselves.

## Writing the bcftools command argument

- DivBase expects bcftools-style command strings
- Only supports `bcftools view`
- Several commands can be piped together by separating them by semicolons
- The backend will apply the commands to each VCF file needed for the query in turn, and finally merge and/concatenate them in to a single results file. This means that the user should not state `merge` or `concat` in their commands.

- samples names are autoinjected when the backend builds the `bcftools` commands. By default `view -s <SAMPLES_FROM_USER_OR_METADATA_QUERY` will be the first command in the pipe. But users can move that the end.

Examples of `bcftools view` subcommands that can be used with DivBase

TODO improve this

|`bcftools view` subcommand, short form | subcommand, long form| Short explanation |
|---|---|--|
|-G | --drop-genotypes| |
|-r | --regions | |
|-A | --trim-unseen-alleles| |
|-a | --trim-alt-alleles| |
|--force-samples| N/A | |
|-I |  --no-update| |
|-s LIST_OF_SAMPLES | --samples LIST_OF_SAMPLES | NOTE! special case in DivBase. Can be used to specify where in a pipe the samples subsetting should occur. Do not specify samples, this is automatically handled by the DivBase server |

TODO: ensure backend strips `-s LIST_OF_SAMPLES` to just `-s`
TODO: since we only support `view`, can there be a shortform where we skip `view` and just have the view flags?

There are more examples of `view` commands that are supported. Not every single one might have been tested. As long as it is not among the blacklisted subcommands below, it can be part of a `--command` string.

Blacklisted `view` subcommands (not supported in DivBase)

| `bcftools view` subcommand | Reason why it is not allowed in DivBase |
|---|---|
|-S, --samples-file FILE | Covered by `divbase-cli query vcf --samples-list`|

Add explicit “invalid patterns” subsection:

Unsupported subcommands names for `view` commands:

- Among the `view` commands, Sample-file flags inside `--command` (`-S` / `--samples-file`) are not allowed. use the CLI `divbase-cli query vcf --samples-file` instead

Examples:

```bash
--command "view -s SAMPLES; view -r 21:15000000-25000000"
```

TODO: now that samples are autoinjected, we need to support no command? a current workaround is to force them to write `view -s`. i.e. how to handle Empty command string

TODO the -s SAMPLES placeholder still lives on in the docs/docstrings

TODO perhaps the default place for `view -s` should be the last place of the command? since it is faster on shorter files?

## What happens after submitting a VCF query? (Job lifecycle and outputs)

### What the user can see after submitting the job

1. CLI returns `Job submitted successfully with task id: <ID>`
2. User monitors status via task-history commands `divbase-cli task-history`
3. On job success, DivBase uploads a result VCF file to project's data storage. If the job fails, user can read the error message in the task-history.
4. User can list/download result files from the project's data storage.

Suggested commands:

```bash
divbase-cli task-history id <TASK_ID>
divbase-cli files ls --project <PROJECT_NAME> --include-results-files
```

### How VCF queries affect files in your project

- Source VCF files: uploaded by project members to the DivBase project. Read, but not modified by the query jobs. Source VCF files are indexed in the VCF dimensions cache of the project.
- Sidecar metadata TSV: read-only during queries.
- Result VCF files: new files created on successful jobs. Are never considered towards the VCF dimensions cache or subsequent VCF queries. Will be named `result_of_job_<JOB_ID>`. As an extra layer of provenence in case the file name of the results file is change, a the following is also added to the file header using `bcftools annotate`: `##DivBase_created="This is a results file created by a DivBase query; Date=<TIME_STAMP>"`

!!! Note
    - Re-running queries creates identical jobs and new result files. The system does not have any limitations for duplicate queries, so please be mindful of this.
    - Result files are outputs, not treated as new source data for future dimensions/indexing workflows
    - Keep bucket tidy by deleting obsolete result files if needed, they will take up space

TODO describe cron job for old results files

### How does DivBase process the VCF files (technical implementation)

subset each input file by itself. apply all commands in the pipe to it. save to a temp file. merge and/or concat all tempfiles into a single results file. upload the results files to the DivBase project.

## Common errors and how to fix them

TODO: Add a short table:

| Error | Likely cause| Suggested fix |
|---|---|---|
|TODO ADD CLI ERROR HERE | TODO ADD CAUSE| TODO ADD FIX |

- “Use only one of --tsv-filter, --samples, or --samples-file”
- Unknown sample IDs
- Empty/invalid `--samples-file` format
- Dimensions out of date
- Invalid filter syntax / unsupported command

See also:

- [Troubleshooting](troubleshooting.md)
- [VCF Dimensions](vcf-dimensions.md)
- [Sample metadata query](sidecar-metadata-queries.md)

## Practical examples

- Region-only query across all samples
- Metadata + VCF combined query
- Direct sample list query
- Multi-step command with semicolon-separated pipeline stages

Each example should include:

- Command
- What it selects
- What output file behavior to expect

## Read next

- [How to create efficient DivBase queries](how-to-create-efficient-divbase-queries.md)
- [Tutorial: Running a query on a public dataset](tutorial-query-on-public-data.md)
