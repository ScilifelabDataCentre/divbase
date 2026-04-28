# e2e tests for VCF query results checksums

In `tests/e2e_integration/cli_commands/test_query_cli.py`, the test `test_vcf_query_result_file_by_headerless_checksum()` is intended to check that the results files of DivBase VCF queries are correct by comparing a checksum of the results file with an expected value. This allows to easily extend the test to other VCF files by adding the query and the expected checksum to the test parameters.

The crux is that VCF files can contain timestamps in the headers (particularly: in the `##` header lines, not in the `#` header line of the VCF file). Timestamps are used for the VCF queries in DivBase: each 'bcftools view' subsetting operation ran to produce the results file during the VCF query workflow is recorded to the `##` header. In many ways, this is a desirable audit trail in DivBase, and whilst it can be turned off, this is not an option for DivBase at the time of writing. Therefore, the parametrised e2e test named `test_vcf_query_result_file_by_headerless_checksum()` in `tests/e2e_integration/cli_commands/test_query_cli.py` has been designed to calculate the `##` headerless checksum for an given vcf query and compare it to its expected `##` headerless checksum.

This document describes how to add new test cases to `test_vcf_query_result_file_by_headerless_checksum()`. The user will need to calculate the expected `##` headerless checksum for a given results VCF file from a query.

## 1. Calculate the expected `##` headerless checksum

The e2e `test_vcf_query_result_file_by_headerless_checksum()` in `tests/e2e_integration/cli_commands/test_query_cli.py` use a helper function `_checksum_vcf_skip_double_hash_headers()` that in turn call a shared function in the codebase for calculating the MD5 checksum from files (`calculate_md5_checksum`). All that is needed to calculate the expected `##` headerless checksum is a results VCF file. The code in [Section 1.1](#11-general-pattern-to-calculate-the--headerless-md5-checksum-for-a-results-vcf-file) can then be used as is from the `divbase` repository root.

Depending on what you want to test, you can either generate the results VCF file with a DivBase VCF query ([Section 1.2](#12-generating-the-results-file-using-divbase-cli-query-vcf)), or by identifying the bcftools workflow sequence for the given query and manually running the command sequence in the Celery worker environment ([Section 1.3](#13-generating-the-results-file-by-manually-running-the-bcftools-steps)).

The latter requires a bit more work and an understanding of how the bcftools orchestration works in DivBase. A middle ground, which is described in [Section 1.3](#13-generating-the-results-file-by-manually-running-the-bcftools-steps) is to run the `divbase-cli query vcf` command, and deduce the bcftools command sequence from the Celery worker logs. A benefit from this method is that it can be a sanity-check to ensure that the bcftools command sequence is created and run as intended. This method is not independent from DivBase since the logs are used the guide the manual run, but it reduces dependence on the layered signal from CLI to API to worker.

### 1.1. General pattern to calculate the ## headerless MD5 checksum for a results VCF file

For a given results VCF file `<PATH_TO_RESULTS_FILE>`, run the following in the terminal to calculated the expected `##` headerless checksum:

```bash
uv run python - <<'PY'
from pathlib import Path
from tempfile import TemporaryDirectory
from tests.e2e_integration.cli_commands.test_query_cli import _checksum_vcf_skip_double_hash_headers

with TemporaryDirectory() as tmp:
    results_vcf_file = "<PATH_TO_RESULTS_FILE>"
    print(_checksum_vcf_skip_double_hash_headers(Path(results_vcf_file), Path(tmp)))
PY
```

### 1.2. Generating the results file using divbase-cli query vcf

Run a VCF query with `divbase-cli query vcf` (assuming that the [prerequisites have been met](../user-guides/vcf-query-syntax.md/#1-prerequisites)). This can be run against a local docker compose stack of DivBase or against a deployed DivBase instance. For the intents of this guide, a local docker compose stack is assumed (`docker compose -f docker/divbase_compose.yaml watch`).

Example (assumes a DivBase project with the fixtures `HOM_20ind_17SNPs_first_10_samples.vcf.gz` and `HOM_20ind_17SNPs_last_10_samples.vcf.gz`):

```bash
divbase-cli query vcf --tsv-filter "Area:West of Ireland,Northern Portugal;Sex:F" --command "view -s"
# Example terminal output
#Job submitted successfully with task id: 281. To check the status of your job, use the command: divbase-cli task-history id 281
```

When the job is finished, download the results file to disk and run the code in [Section 1.1](#11-general-pattern-to-calculate-the--headerless-md5-checksum-for-a-results-vcf-file) in your terminal.

For the above example: if the results file was downloaded to the root of the `divbase` repo, `<PATH_TO_RESULTS_FILE>` would be `result_of_job_281.vcf.gz`.

The checksum for the example should be:

```bash
3f9c371bcffb8126663cf08a802ae58c
```

### 1.3. Generating the results file by manually running the bcftools steps

An alternative to generating the results file from DivBase itself is to manually recreate the bcftools steps generated by the DivBase backend when processing the query. This can quickly become complex since DivBase can act on multiple source VCF files and run several subsequent commands on each file. To see exactly how DibBase processes a given VCF query, you can submit a job and inspect the Celery worker logs after the job has finished.

Example:

```bash
divbase-cli query vcf --tsv-filter "Area:West of Ireland,Northern Portugal;Sex:F" --command "view -s"
```

Check the Celery worker logs. VCF queries are wired to be processed by the `divbase-worker-long` container in the local Docker Compose stack of DivBase.

```bash
# Assuming the you ran the query against a local docker compose stack of DivBase
docker logs divbase-worker-long-1
```

From the logs, identify the `bcftools` commands, their order and their options. Take a moment to analyse if these are reasonable commands. This is a good sanity-check that can help find bugs in the orchestration system.

It is desirable to run the manual `bcftools` in a DivBase worker container to ensure that versions of all tools and dependencies. Later on when the VCF case is to be implement as a parameter for the e2e test, the VCF files need to be stored under source control in `tests/fixtures`. The testing docker compose stack already is configured to mount this directory, so we can use that to get easy access to the source VCF files needed for the query (in the below example: `HOM_20ind_17SNPs_first_10_samples.vcf.gz` and `HOM_20ind_17SNPs_last_10_samples.vcf.gz`).

```bash
# Spin up the testing stack since it mounts the tests/fixture dir as per docker/divbase_compose.tests.yaml
docker compose -f docker/divbase_compose.yaml -f docker/divbase_compose.tests.yaml up -d --build

docker exec -it -w /app/tests/fixtures divbase-tests-worker-long-1 sh

# Inside container, run the bcftools commands:
# For the example 'divbase-cli query vcf --tsv-filter "Area:West of Ireland,Northern Portugal;Sex:F" --command "view -s"'
bcftools view -s 5a_HOM-I7,1b_HOM-G58 HOM_20ind_17SNPs_first_10_samples.vcf.gz -Ou -o temp_subset_0_0.bcf
bcftools index -f temp_subset_0_0.bcf

bcftools view -s 5a_HOM-I13,5a_HOM-I14,5a_HOM-I20,5a_HOM-I21 HOM_20ind_17SNPs_last_10_samples.vcf.gz -Ou -o temp_subset_0_1.bcf
bcftools index -f temp_subset_0_1.bcf

bcftools merge --force-samples -Ou -o merged_unsorted.bcf temp_subset_0_0.bcf temp_subset_0_1.bcf

printf '##DivBase_created="This is a results file created by a DivBase query; Date=Mon Apr 27 09:37:12 2026"\n' > divbase_header.txt
bcftools annotate -h divbase_header.txt -Ou -o merged_annotated_unsorted.bcf merged_unsorted.bcf
bcftools sort -Oz -o result_of_job.vcf.gz merged_annotated_unsorted.bcf

exit

# Back on host:
uv run python - <<'PY'
from pathlib import Path
from tempfile import TemporaryDirectory
from tests.e2e_integration.cli_commands.test_query_cli import _checksum_vcf_stream_skip_double_hash_headers

with TemporaryDirectory() as tmp:
 results_vcf_file="tests/fixtures/result_of_job.vcf.gz"
 print(_checksum_vcf_stream_skip_double_hash_headers(Path(results_vcf_file).read_bytes(), Path(tmp)))
PY

# Spin down the testing stack
docker compose -f docker/divbase_compose.yaml -f docker/divbase_compose.tests.yaml down -v
```

For this example, this should result in the checksum:

```bash
3f9c371bcffb8126663cf08a802ae58c
```

## 2. Implementing the test case in the parametrized test for ## headerless VCF checksums
