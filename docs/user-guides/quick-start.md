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

Managers can see the guide on [Account management](account-management.md) for details on how to add members to a project.

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

## Step 4: Add your project(s) to your divbase-cli config

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
divbase-cli files ls
```

## Step 7: Dimensions update

For DivBase to be able to efficiently handle the VCF files in the the project, some key information about each VCF files is fetched from the files and cached in the server. In DivBase, this is refered to as "VCF dimensions". These include for instance which samples and scaffolds that a VCF file contains.

Update the project dimensions after uploading your files.

```bash
divbase-cli dimensions update
```

This submits a task to the DivBase task management system. The task will wait in a queue until the system is ready to work on it. Depending on the size of the VCF files, this make take a couple of minutes.

!!! notes
    1. Please note that it is not possible to run VCF queries in DivBase until the dimensions update task has finished. The reason for this is that the VCF queries use the dimensions data ensure that the queries are feasible and to know which VCF files from the project to process.

    2. Please also note that the `divbase-cli dimensions update` command needs to be run every time a new VCF or a new version of a VCF file is uploaded.

## Step 8: Confirm dimensions update job completion

Check the task history to confirm the dimensions update job has completed:

```bash
divbase-cli task-history user
```

Once complete, you can run any queries on the uploaded data.

It is possible to inspect the cached VCF dimensions data for the project at any time with the command:

```bash
divbase-cli dimensions show
```

## Step 9: Upload sample metadata

DivBase can checkout data based the VCF files themselves, but can also take an optional sidecar sample metadata file into account. The metadata file must be a TSV (tab-separated variables) file. The metadata contents of the file is defined by the users. If the VCF dimensions command has been run for the project, the cached dimensions data can be used create a template where the samples of the project have been pre-filled:

```bash
divbase-cli dimensions create-metadata-template
```

Details on how to write this file are given in [Sidecar Metadata TSV files: creating and querying sample metadata files](sidecar-metadata.md). In short, the first row starts with `#` and contains the headers for different metadata columns. The first column (`Sample_ID`) is mandatory and can be created by the system as just described; if created manually just make sure that each sample name is spelled exactly as in the VCF files. The rest of the columns are free for the user to define.

Example of a sidecar metadata TSV file with the mandatory `Sample_ID` column and two user defined columns.

```
#Sample_ID Population Area
129P2 1 North
129S1 2 East
129S5 3 South
```

!!! note
    Please use a text editor that preserves the tabs when the file is saved. Incorrect tabs can lead to issues with running metadata queries in DivBase.

There is a command to help check that the sidecar metadata TSV is correctly formatted for use with DivBase. Running it is optional:

```bash
divbase-cli dimensions validate-metadata-file path/to/your/sample_metadata.tsv
```

When you are happy with the sample metadata file, it should be uploaded to the DivBase project with the following:

```bash
divbase-cli files upload path/to/your/sample_metadata.tsv
```

## Step 10: Run your queries

There are three types of queries in DivBase:

- Sample metadata query
- VCF data query
- Combined sample metadata and VCF data query

!!! note
    Queries are one of the more complex aspects of DivBase and therefore the user is encouraged to read the section on [Running Queries](running-queries-overview.md) after reading this quick start.

### Running sample metadata queries

As an example, let's assume that user-defined sidecar sample metadata TSV file contains a custom column name `Area` and that `Northern Portugal` is one of the values in the column. To filter on all samples on that column and value:

```bash
divbase-cli query tsv "Area:Northern Portugal"
```

!!! note
    Please see [Sidecar Metadata TSV files: creating and querying sample metadata files](sidecar-metadata.md) for more details on the syntax for writing sample metadata queries.

### Running VCF data queries

DivBase uses [`bcftools`](https://github.com/samtools/bcftools) to subset VCF data and therefore the query syntax is based on `bcftools view` syntax (as described in the [bcftools manual](https://samtools.github.io/bcftools/bcftools.html#view)).

For instance, to subset all VCF files in the project on a chromosomal region in a scaffold named `21`:

```bash
divbase-cli query vcf --command "view -r 21:15000000-25000000"
```

The VCF queries can be combined with sidecar sample metadata queries using `--tsv-filter`. DivBase will first resolve which samples match the metadata filter, then automatically inject those sample IDs into the `bcftools` command. Only VCF files that contain at least one of the matching samples will be processed. An example:

```bash
divbase-cli query vcf --tsv-filter "Area:Northern Portugal" --command "view -r 21:15000000-25000000"
```

By default, sample IDs are injected at the first command segment (e.g. the above becomes `view -s S1,S2 -r 21:15000000-25000000`). If you are using a multi-segment pipe and want explicit control over which segment receives the sample injection, you can use `view -s` (without sample values) as a placeholder in that position:

```bash
divbase-cli query vcf --tsv-filter "Area:Northern Portugal" --command "view -s; view -r 21:15000000-25000000"
```

!!! note
    DivBase only allows `bcftools view` in its query syntax and no other `bcftools` commands. The `merge`, `concat`, and `annotate` commands are used when processing a query, but should not be defined by the user.

## Step 11: Download any results files

You can check the status of all of your submitted jobs using:

```bash
divbase-cli task-history user
```

Once a `vcf` job is complete, you can download the resulting merged vcf file:

```bash
divbase-cli files download result_of_job_<JOB_ID>.vcf.gz # --download-dir path/to/save/results/
```

Replacing <JOB_ID> with the actual job ID from the task history.

## Next steps

This quick start guide was hopefully enough to get you started using DivBase. If you want to learn more about the different features of DivBase, we have many detailed guides.

For details on:

- Installing the DivBase client on your computer, see [Installation](installation.md) and [Setup DivBase CLI](setup_divbase_cli.md)

- How to format your VCF files to get the most out of DivBase, see [Working with VCF Files in DivBase](vcf-files.md)

- Everything releated to DivBase queries, we reccomend to start at [Running Queries: Overview](running-queries-overview.md)

- Creating a snapshot the version of the files in a DivBase project at a current time, see [Project versioning](project-versioning.md)

## Getting help

If you run into issues:

1. **Check the output**: Read the error message
2. **Consult and search the docs**: Browse the [Full documentation](../index.md) or look at common issues in the [Troubleshooting](troubleshooting.md) page.
3. **Get help on the command you're running**: `divbase-cli COMMAND --help`

To get assistance from us you can either send us an email (<TODO@scilifelab.se>) or report an issue on our [GitHub Issues](https://github.com/ScilifelabDataCentre/divbase/issues).
