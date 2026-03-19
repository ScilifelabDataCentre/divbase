# DivBase VCF query syntax

Users can checkout subsets of their VCF data from their DivBase project using the command `divbase-cli query vcf`. This data checkout is run as an asynchronous job that is sent to the queue on the DivBase server, and eventually run once there are idle resources to process the job. Users will get a job ID when they submit their query to the task queue, and can view the status of the job with for instance `divbase-cli task-history id <JOB_ID>`.

The processing of the VCF files on the DivBase server is done with [`bcftools`](https://github.com/samtools/bcftools). DivBase will detect the VCF files in the project's data store that are needed for the query; if more than one VCF file is needed, DivBase will ensure that the files are compatible with each other according to the requirements of `bcftools` and ensure that a single results file with the subset data is returned to the user by running `bcftools merge` and `bcftools concat` as needed. The result is a single VCF file that is uploaded to the projects data store and named after the job ID.

Users can query the VCF data in their project with or without combining it to a [sample metadata query](docs/user-guides/sidecar-metadata.md).

Example of a VCF query that identifies the samples and VCF files to filter on in the project's datastore and then applies a subset based on genomic range:

```bash
divbase-cli query vcf --tsv-filter  "Area:North,West;Weight:>10"" --command "view -r 21:15000000-25000000"
```

The outcome a DivBase VCF query is a single results file with merged/concatenated data that fulfills all filters.

!!! Note
    When you submit a query, DivBase will use the state of the VCF Dimensions and the VCF files at that very point in time to produce the query results. It is therefore fine if you or another project member uploads new VCF files to the project while a query is queued or running.

## 1. Prerequisites

DivBase uses [`bcftools`](https://github.com/samtools/bcftools) to subset VCF data and therefore the query syntax is based on `bcftools view` syntax (as described in the [bcftools manual](https://samtools.github.io/bcftools/bcftools.html#view)).

bcftools rules

## 2. `divbase-cli query vcf` command structure

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

## 3. Sample-selection

### 3.1. Mode A: Metadata-driven selection (--tsv-filter)

- Combined sample metadata and VCF data query

User passes a sample metadata query that identifies the samples and VCF files that are needed for the VCF query.

```bash
divbase-cli query vcf --tsv-filter  "Area:North,West;Weight:>10"" --command "view -r 21:15000000-25000000"
```

For example, the system might find that the samples that fulfill the [sample metadata query](docs/user-guides/sidecar-metadata.md) set with `--tsv-filter` are, say, `S2`, `S5` `S28` `S108` and that they described in the files `file1.vcf`, `file3.vcf`, `file4.vcf`. The DivBase server will then act on only these three files and subset based on the four samples.

### 3.2. Mode B1: Direct sample list (--samples)

- Regular VCF data query (without metadata query)

To get all samples that are available for a project, use

User specifies the samples to query.

```bash
divbase-cli query vcf --samples "S1,S2,S10,S239" --command "view -r 21:15000000-25000000"
```

it IS possible to specify both the samples and the files

```bash
divbase-cli query vcf --samples "S1,S2,S10,S239" --vcf-files "file1.vcf, file5.vcf" --command "view -r 21:15000000-25000000"
```

but this requires the user to know which samples are in which file. If

--vcf-files "file1.vcf, file5.vcf"

### 3.3. Mode B2: Samples from file (--samples-file)

```bash
divbase-cli query vcf --samples-file samples_for_my_query.txt --command "view -r 21:15000000-25000000"
```

The file needs to be a plain file, such a `.txt` and contain one sample per row.

```text
S1
S2
S10
S239
```

### 3.4. Mode C: No selection flags (all samples)

TODO: we might not want to support this, it would be easier for the user to download all the files and merge them themselves.

## 4. Writing the `--command` argument

- DivBase expects bcftools-style command strings
- Only supports `bcftools view`
- Several commands can be piped together by separating them by semicolons
- The backend will apply the commands to each VCF file needed for the query in turn, and finally merge and/concatenate them in to a single results file. This means that the user should not state `merge` or `concat` in their commands.

Add explicit “invalid patterns” subsection:

Unsupported subcommands names for `view` commands:

- Among the `view` commands, Sample-file flags inside `--command` (`-S` / `--samples-file`) are not allowed. use the CLI `divbase-cli query vcf --samples-file` instead

Examples:

```bash
--command "view -s SAMPLES; view -r 21:15000000-25000000"
```

TODO: now that samples are autoinjected, we need to support no command? a current workaround is to force them to write `view -s`. i.e. how to handle Empty command string

TODO the -s SAMPLES placeholder still lives on in the docs/docstrings

## 5. What happens after submitting a VCF query? (Job lifecycle and outputs)

### 5.1. What the user can see after submitting the job

1. CLI returns `Job submitted successfully with task id: <ID>`
2. User monitors status via task-history commands
3. On success, DivBase writes a result VCF file to project storage
4. User can list/download result files

Suggested commands:

```bash
divbase-cli task-history id <TASK_ID>
divbase-cli files ls --project <PROJECT_NAME> --include-results-files
```

### 5.2. How VCF queries affect files in your project

Distinguish between:

- Source VCF files: not modified by query jobs
- Sidecar metadata TSV: read-only during query
- Result VCF files: new files created on successful jobs

Note!

- Re-running queries creates additional result files
- Result files are outputs, not treated as new source data for future dimensions/indexing workflows
- Keep bucket tidy by deleting obsolete result files if needed, they will take up space

TODO describe cron job for old results files

### 5.3. How does DivBase process the VCF files (technical implementation)

subset each input file by itself. apply all commands in the pipe to it. save to a temp file. merge and/or concat all tempfiles into a single results file. upload the results files to the DivBase project.

## 6. Common errors and how to fix them

TODO: Add a short table:

- Symptom
- Likely cause
- Fix

Rows

- “Use only one of --tsv-filter, --samples, or --samples-file”
- Unknown sample IDs
- Empty/invalid `--samples-file` format
- Dimensions out of date
- Invalid filter syntax / unsupported command

Cross-link to:

- `troubleshooting.md`
- `vcf-dimensions.md`
- `sidecar-metadata-queries.md`

## 8. Practical examples

- Region-only query across all samples
- Metadata + VCF combined query
- Direct sample list query
- Multi-step command with semicolon-separated pipeline stages

Each example should include:

- Command
- What it selects
- What output file behavior to expect
