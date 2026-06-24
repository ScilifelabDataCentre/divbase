# Tips for using DivBase programmatically

Below is a set of tips for users who want to use DivBase in scripts/pipelines/programmatically.

If you have any tips or suggestions to add to this page or any desired features, please let us know!

## Use Personal Access Tokens to authenticate programmatically

For scripts, pipelines, and HPC jobs the recommended approach is to use a **Personal Access Token (PAT)**. A PAT is a static bearer token that you create once via the website. When using a PAT, there is no login step and no password storage required.

See [Account Management — Personal Access Tokens](./account-management.md#personal-access-tokens) for how to create and revoke PATs.

You can store/use a PAT in two ways

### (Recommended) — Let `divbase-cli` handle storing the PAT

After creating a PAT on the DivBase website, store it on your device by copy pasting the pre-filled commands shown on the website, it will look something like this:

```bash
# with an expiry date:
divbase-cli auth add-pat "my-pat-name" --expires UNIX_TIMESTAMP
# or without an expiry date:
divbase-cli auth add-pat "my-pat-name"
```

You will be prompted to paste the token value. It will then be stored securely in your OS keyring (or in a restricted file if no keyring is available) and used automatically on every subsequent `divbase-cli` command.

```bash
divbase-cli auth pat-info   # show name and expiry of the stored PAT
divbase-cli auth rm-pat     # remove the stored PAT from this device
```

!!! info "This strategy works on HPC clusters"
    On the login node with divbase-cli [installed as we recommend](./installation.md), store the PAT as described above. The token will be available to all your jobs running on the cluster, without needing to worry about setting any environment variables in your job scripts.

### (Not recommended alternative) — Store the PAT in an environment variable

Set `DIVBASE_API_PAT` in your shell or job script:

```bash
export DIVBASE_API_PAT="divbase_pat_your_token_here"
divbase-cli files ls
```

You may prefer this if you want explicit per-job control of which token is used.

---

!!! question "What if I have both an active login session and a personal access token set?"
    `divbase-cli` follows this priority order:

    1. **Active login session** (from `divbase-cli auth login`) — highest priority
    2. **`DIVBASE_API_PAT` environment variable**
    3. **CLI-stored PAT** (from `divbase-cli auth add-pat`)

    To use a PAT when you are also logged in, run `divbase-cli auth logout` first.

    Your logged in account will have the same or more permissions that your PATs, this is why we follow this order.

## Parse divbase-cli files ls/info output programmatically

You can make the output of the `divbase-cli files info` and `divbase-cli files ls --detailed` commands in TSV format for easier parsing. Use the `--tsv` flag:

    ```bash
    divbase-cli files ls --tsv
    divbase-cli files info FILE_NAME --tsv
    ```

    You can do the same for any [project versions](./project-versioning.md) you've created for your project:

    ```bash
    divbase-cli version ls --tsv
    divbase-cli version info VERSION_NAME --tsv
    ```

## Stream files from the command line

Rather than first downloading a file, you can stream a file from the command line and pipe it into other tools for processing directly without saving it to disk.

    ```bash
    divbase-cli files stream my_file.vcf.gz | zcat | less
    ```

!!! info "BCFTools accepts stdin as input"
    You can pipe a VCF file directly into BCFTools without downloading it first:

    ```bash
    divbase-cli files stream my_file.vcf.gz | bcftools view -h -
    ```

    Especially useful if you just want to view the header of a VCF file without having to download the entire file.

## Running VCF queries programmatically

VCF queries can potentially take some time depending on the length of the queue on the DivBase server and the time required to process the query. The `divbase-cli query get-vcf-results` CLI command has been implemented with scripts and automated worksflows in mind: it polls the submitted task for completion, and if the task finishes with `SUCCESS`, it downloads the results file for further processing.

Example bash script for submitting a VCF query with desired `<parameters>`, polling for its successful completion and, if successful, downloading of the results file to `<path/to/destination/directory>`:

```bash

# Submit a VCF query with the desired <parameters> and extract the task ID from the terminal output upon successful submission
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
    # task ID does not belong to a VCF query task
fi
```

!!! Note
    It is also possible to check for the status of a submitted task with `divbase-cli task-history id <task_id>`, but that command does not have built in polling (and downloading) like `divbase-cli query get-vcf-results` does.
