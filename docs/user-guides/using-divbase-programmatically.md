# Tips for using DivBase programmatically

Below is a set of tips for users who want to use DivBase in scripts/pipelines/programmatically.

If you have any tips or suggestions to add to this page or any desired features, please let us know!

## Use Personal Access Tokens to Authenticate programmatically

For scripts, pipelines, and HPC jobs the recommended approach is to use a **Personal Access Token (PAT)**. A PAT is a static bearer token that you create once via the website and pass to `divbase-cli` via an environment variable. When using a PAT, there is no login step and no password storage required.

See [Account Management — Personal Access Tokens](./account-management.md#personal-access-tokens) for how to create/and remove PATs.

Once you have a token, set the `DIVBASE_API_PAT` environment variable to it. `divbase-cli` will automatically use it in every request.

!!! question "What if I have both an active login session and a Personal Access Token set?"
    `divbase-cli` prioritises an active login session over a PAT. If you have both, the CLI will use the active session and ignore the PAT. To use the PAT, you would need to run `divbase-cli auth logout` first.

```bash
export DIVBASE_API_PAT="divbase_pat_your_token_here"
divbase-cli files ls
```

When `DIVBASE_API_PAT` is set, `divbase-cli` does not need you to be logged in.

### Example: Slurm job script

The cleanest way to use a PAT in a SLURM job is to store the token in a restricted file and load it at job start:

```bash
echo "divbase_pat_your_token_here" > ~/.divbase_pat
chmod 600 ~/.divbase_pat  # only readable/writeable by the owner
```

Then in your SLURM script:

```bash
#!/bin/bash
#SBATCH --job-name=my_divbase_job
#SBATCH --time=24:00:00
# ....

export DIVBASE_API_PAT=$(cat ~/.divbase_pat)

# Download the files you need
divbase-cli files download my_data.vcf.gz
```

!!! tip "Scope your token to what the job needs and when you need it for"
    When creating the PAT, restrict it to the specific project(s) you need it for. Consider also setting an appropriate expiry date for the token. You can always revoke the token immediately if needed from the DivBase website.

## Parse divbase-cli files ls/info output programmatically

1. You can make the output of the `divbase-cli files info` and `divbase-cli files ls` commands in TSV format for easier parsing. Use the `--tsv` flag:

    ```bash
    divbase-cli files ls --tsv
    divbase-cli files info FILE_NAME --tsv
    ```

    You can do the same for any [project versions](./project-versioning.md) you've created for your project:

    ```bash
    divbase-cli version ls --tsv
    divbase-cli version info VERSION_NAME --tsv
    ```

2. Rather than first downloading a file, you can stream a file from the command line and pipe it into other tools for processing directly without saving it to disk.

    ```bash
    divbase-cli files stream my_file.vcf.gz | zcat | less
    ```

    !!! Info
        BCFTools accepts stdin as input, so you can also pipe a VCF file directly into BCFTools without saving it first:

        ```bash
        divbase-cli files stream my_file.vcf.gz | bcftools view -h -
        ```

## Running VCF queries programmatically

VCF queries can potentially take some time depending on the length of the queue on the DivBase server and the time required to process the query. The `divbase-cli query get-vcf-results` CLI command has been implemented with scripts and automated worksflows in mind: it polls the submitted task for completion, and if the task finishes with `SUCCESS`, it downloads the results file for further processing.

Example bash script for submitting a VCF query with desired `<parameters>`, polling for its successful completion and, if successful, downloading of the results file to `<path/to/destination/directory>`:

```bash

# Submit a VCF query with the desired <parameters> and extract the DivBase Task ID from the terminal output upon successful submission
TASK_ID=$(divbase-cli query vcf <parameters> | sed -n 's/.*task id: \([0-9]*\)\..*/\1/p')
if [ "${PIPESTATUS[0]}" -ne 0 ]; then
    echo "Error: failed to submit VCF query job. See output above for details." >&2
    exit 1
fi

divbase-cli query get-vcf-results "$TASK_ID" --download-dir <path/to/destination/directory>
EXIT_CODE=$?
if [ $EXIT_CODE -eq 0 ]; then
    # results file downloaded to <path/to/destination/directory>, process it
elif [ $EXIT_CODE -eq 1 ]; then
    # task failed with Celery task status FAILURE
elif [ $EXIT_CODE -eq 2 ]; then
    # Task ID does not belong to a VCF query task
fi
```

!!! Note
    It is also possible to check for the status of a submitted task with `divbase-cli task-history id <TASK_ID>`, but that command does not have built-in polling (and downloading) like `divbase-cli query get-vcf-results` does.
