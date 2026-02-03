# Quick Start Guide

This guide will walk you through everything you need to start managing and querying your VCF files.

## Prerequisites

- Python 3.12 or higher
- VCF files and sample metadata in TSV format

## Step 1: Create an account

Create an account on [DivBase](https://divbase.scilifelab.se) and make sure to verify your email address.

## Step 2: Join or create a project

### Option 1: Join an existing project

If you have an account you can be added to an existing project by a project member with the role manager. Ask them to add you and give them the email address you used to sign up with.

Managers can see this guide here on how to add members to a project here (TODO).

### Option 2: Create a new project

To create a new project, you'll need to contact us at <TODO@scilifelab.se> with the following information:

- Desired project name
- A brief description of the project
- Your registered email address

## Step 3: Install divbase-cli

Install `divbase-cli` using pipx (recommended):

```bash
pipx install divbase-cli
```

If you do not have `pipx` installed, you can install it by following [the official instructions from pipx](https://pipx.pypa.io/stable/installation/). Refer to the [Installation Guide](installation.md) for more detailed instructions or other ways to install divbase-cli.

## Step 4: Configure the CLI

Set up your user configuration file:

```bash
divbase-cli config create
```

This creates a configuration file stored on your local device at: `~/.config/divbase/config.yaml`.

Add your project(s) to the configuration file:

```bash
divbase-cli config add PROJECT_NAME --default
```

By setting this as the default project, you won't need to specify the project name in future commands. To select a different project in any future command, you can add the `--project` flag to any command. Note that the rest of the example commands in this quick start guide will omit the `--project` flag for brevity.

!!! note
    On the divbase website you can see your project's name(s) under the "Projects" tab after logging in.

## Step 5: Log in

Log in to DivBase:

```bash
divbase-cli auth login EMAIL_ADDRESS
```

## Step 6: Upload files

Upload your VCF files to your project:

```bash
# Upload a single VCF file
divbase-cli files upload path/to/your/file.vcf

# Upload multiple files
divbase-cli files upload path/to/file1.vcf path/to/file2.vcf

# Upload all VCF files in a directory
divbase-cli files upload --upload-dir /path/to/directory/
```

Check your uploaded files:

```bash
divbase-cli files list
```

## Step 7: Upload sample metadata

Sample metadata must be uploaded as follows:

- In TSV format and be named "sample_metadata.tsv"
- Must contain a column named "sample_id" which matches the sample IDs in your VCF files
- The names and values of all other columns are optional.

TODO update this after `sidecar-metadata.md` docs are done, there are changes planned for some details.

```bash
divbase-cli files upload path/to/your/sample_metadata.tsv
```

## Step 8: Dimensions update

For DivBase to be able to efficiently handle the VCF files in the the project, some key information about each VCF files is fetched from the files. In DivBase, this is refered to as "VCF dimensions". These include for instance which samples and scaffolds that a VCF file contains.

Update the project dimensions after uploading your files.

```bash
divbase-cli dimensions update
```

This submits a task to the DivBase task management system. The task will wait in a queue until the system is ready to work on it. Depending on the size of the VCF files, this make take a couple of minutes.

!!! notes
    1. Please note that it is not possible to run VCF queries in DivBase until the dimensions update task has finished. The reason for this is that the VCF queries use the dimensions data ensure that the queries are feasible and to know which VCF files from the project to process.

    2. Please also note that the `divbase-cli dimensions update` command needs to be done every time a new VCF or a new version of a VCF file is uploaded.

## Step 9: Confirm dimensions update job completion

Check the task history to confirm the dimensions update job has completed:

```bash
divbase-cli task-history user
```

Once complete, you can run any queries on the uploaded data.

It is possible to inspect the cached VCF dimensions data for the project at any time with the command:

```bash
divbase-cli dimensions show
```

## Step 10: Run your queries

There are three types of queries in DivBase:

- Sample metadata query
- VCF data query
- Combined sample metadata and VCF data query

!!! notes
    Queries are one of the more complex aspects of DivBase and therefore the user is encouraged to read the section on [Running Queries](running-queries.md) after reading this quick start.

### Running sample metadata queries

As an example, let's assume that user-defined sidecar sample metadata TSV file contains a custom column name `Area` and that `Northern Portugal` is one of the values in the column. To filter on all samples on that column and value:

```bash
divbase-cli query tsv "Area:Northern Portugal"
```

!!! notes
    Please see [Sidecar Metadata TSV files: creating and querying sample metadata files](sidecar-metadata.md) for more details on the syntax for writing sample metadata queries.

### Running VCF data queries

DivBase uses [`bcftools`](https://github.com/samtools/bcftools) to subset VCF data and therefore the query syntax is based on `bcftools view` syntax (as described in the [bcftools manual](https://samtools.github.io/bcftools/bcftools.html#view)).

For instance, to subset all VCF files in the project on a chromosomal region in a scaffold named `21`:

```bash
divbase-cli query bcftools-pipe --command "view -r 21:15000000-25000000"
```

The VCF queries can be combined with sidecar sample metadata queries with `--tsv-filter` and the fixed expression `view -s SAMPLES` (where `SAMPLES` tells DivBase to use the results from the sidecar filtering as input for `bcftools view -s`). In this way, only the VCF files that fulfil the sample metadata query will be used in the `bcftools` subset commands. An example:

```bash
divbase-cli query bcftools-pipe --tsv-filter "Area:Northern Portugal" --command "view -s SAMPLES; view -r 21:15000000-25000000"
```

## Step 11: Download any results files

You can check the status of all of your submitted jobs using:

```bash
divbase-cli task-history user
```

Once a `bcftools-pipe` job is complete, you can download the resulting merged vcf file:

```bash
divbase-cli files download merged_[JOB_ID].vcf.gz # --download-dir path/to/save/results/
```

Replacing [JOB_ID] with the actual job ID from the task history.

## Next steps

TODO - a selection of links to more detailed user guides

## Getting help

If you run into issues:

1. **Check the output**: Read the error message
2. **Consult and search the docs**: [Full documentation](../index.md)
3. **Get help on the command you're running**: `divbase-cli COMMAND --help`

To get assistance from us you can either send us an email (TODO - link) or report an issue on our [GitHub Issues](https://github.com/ScilifelabDataCentre/divbase/issues).

## Common Issues

??? question "Authentication Issues"
    - Make sure you've verified your email address
    - Check if you can login to your account on the [DivBase Website](https://divbase.scilifelab.se). If it fails on the website it will also fail on the CLI.
    - Try logging out and back in: `divbase-cli auth logout` then `divbase-cli auth login`
