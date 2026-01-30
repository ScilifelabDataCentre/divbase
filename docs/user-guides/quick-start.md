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

By setting this as the default project, you won't need to specify the project name in future commands. To select a different project in any future command, you can add the `--project` flag to any command.

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

TODO - do we support numerical columns, including queries on them, i.e in range 31 - 50. etc...

```bash
divbase-cli files upload path/to/your/sample_metadata.tsv
```

## Step 8: Dimensions update

Update the project dimensions after uploading your files:

```bash
divbase-cli dimensions update
```

## Step 9: Confirm dimensions update job completion

Check the task history to confirm the dimensions update job has completed:

```bash
divbase-cli task-history user
```

Once complete, you can run any queries on the uploaded data.

## Step 10: Run your queries

Query your VCF data:

```bash
# TODO - some simple examples to show what kind of queries can be run
```

!!! note
    For more detailed information on running queries, refer to the [Running Queries Guide](running-queries.md).

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
