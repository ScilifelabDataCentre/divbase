# DivBase

**DivBase's goal is to be a service which is comprised of:**

- An S3 object store to store researchers VCF files and associated metadata.
- A CLI installed by users to interact with their DivBase projects and the webAPI (next bullet).
- A webAPI, installed on SciLifeLab cluster to run query jobs and create "results files" from these queries.
  - The results files are added to the object store by and can be downloaded by the researchers.
- As you may have guessed, the CLI can therefore requests to the webAPI including making queries.submitting query jobs and can download/upload files from the object store.

## Table of Contents

1. [Folders overview](#folders-overview)
2. [Deployment of DivBase](#deployment-of-divbase)
3. [divbase CLI](#divbase-cli)
4. [Developer Setup](#developer-setup)
5. [Queries](#queries)

## Folders overview

- *docker:* bcftools docker image.
- *docs:* Documentation
- *src:* Source code for divbase's CLI and webAPI.
- *tests*: Tests (with Pytest).

## Deployment of DivBase

DivBase's deployment with k8s is managed in our [private repository, argocd-divbase](https://github.com/ScilifelabDataCentre/argocd-divbase)

## divBase CLI

The CLI tool `divbase-cli` is designed to help you interact with your DivBase project. The CLI is built using [Typer](https://typer.tiangolo.com).

- See the basic usage/documentation at [docs/using_divbase_cli.md](docs/using_divbase_cli.md)
- See the [Developer Setup](#developer-setup) section to install the CLI tool in develop mode.

## Developer Setup

### 1. Setup your python virtual environment

You can use either [uv](https://github.com/astral-sh/uv) or something like `venv` or `pip`.

#### Using `uv` (Recommended)

- Install `uv`, follow the docs: <https://docs.astral.sh/uv/>

- In the root of the repository, run:

```bash
uv sync
```

This will create a virtual environment at the root of the repo and install all dependencies and the package (with all code in the src folder installed in "editable" mode).

#### Using pip and venv

```bash
python3 -m venv .venv
source .venv/bin/activate
# The `-e .` flag installs the package in "editable" mode.
pip install -e .
```

#### Available commands

You can now run the cli tool and webAPI/server with the following commands.

```bash
divbase-cli
divbase-api # NOTE: We typically do not run the API like this, instead we use docker compose, keep reading below to see how.
```

### 2. Install pre-commit hooks

We also use [pre-commit hooks](https://pre-commit.com/). pre-commit runs on every commit.

```bash
pre-commit install
```

### 3. CLI Setup

divbase-cli relies on 2 local files to preserve state between runs:

1. A config file stored in your home directory at `~/.config/divbase/config.yaml`.
    This config contains information about your DivBase projects and which divbase server (if any) you're logged into.
2. A secrets file stored in your home directory at `~/.config/divbase/.secrets`.
    This file contains your access/refresh tokens for the DivBase server you're logged into.

You can read the [docs/using_divbase_cli.md](docs/using_divbase_cli.md) for more information on how to setup the CLI tool for first time use.

For convienance, we provide a local development setup script, but in order to run this you need the docker compose stack running (see next section).

### 4. Run DivBase server locally using Docker compose watch

```bash
docker compose -f docker/divbase_compose.yaml down && docker compose -f docker/divbase_compose.yaml watch
```

This will deploy the API, job system and MinIO (S3) instance locally. Using `compose watch` means changes in the `/src` folder will automatically be sycned to the container and if needed trigger a restart or rebuild of them.

Once the backend is up and running you can then run `divbase-cli` commands against it or go to for example <http://localhost:8000/api/v1/docs> to test out some of API commands directly.

There is now a very simple DivBase frontend running on <http://localhost:8000>. Frontend is part of the FastAPI app deployment. All API routes are appended with /api/v1 to avoid collisions. Frontend returns HTML etc... directly to browser via Jinja2 templating.

To help with setup, we provide a local development setup script which will:

- Create buckets on MinIO and add some data.
- Create users and projects on DivBase and add users with different roles to each project.
- create a user config file with these projects added to it.

You can run the local_dev_setup script like this:

```bash
uv run scripts/local_dev_setup.py # make sure the compose stack is running
```

### 5. Running tests

- Tests are written using Pytest. Tests are split at the top folder level into `unit` and `e2e_integration` tests (alongside some shared fixtures).
- The distinction between `unit` and `e2e_integration` here is any test that requires some/all of the docker compose testing stack to be running is considered an end-to-end/integration test.
- Note that there maybe "unit like" tests inside the `e2e_integration` dir, primarily for convenience/legacy reasons (before this top level split was made).
- Unit tests should follow the pattern `tests/unit/<package>/<module>/<test_file_name>.py` (i.e. match the module structure of the packages in the codebase).

To run all tests:

```bash
pytest # you may want to append the -s flag to print standard output.
```

To run only unit or e2e_integration tests:

```bash
pytest tests/unit
pytest tests/e2e_integration
```

**The e2e_integration tests will be slower the first time you run them as the docker images will need to be downloaded and built. If you use "-s" you'll see the status of the docker compose building steps.**

## Queries

The data and metadata stored in a DivBase project can be queried two main ways: by a user-submitted sample metadata file, and by the VCF files. The examples in this section intend to show how a few of the commands in `divbase-cli query` can be used. Please add the `--help` flag to any of the commands for more documentation on their usage and options.

### 1. Before performing any queries: update the VCF dimensions file

To improve the speed and robustness of the queries in DivBase, technical metadata from the VCF files in the buckets are stored in the backend in a so-called "VCF dimensions" file. At the moment, users need to run the following command to ensure that the dimensions file is up-to-date with the current VCF files contained in the S3 Object Store of the project:

```bash
divbase-cli dimensions update --project <PROJECT_NAME>
```

To display the current status of the dimensions file:

```bash
divbase-cli dimensions show --project <PROJECT_NAME>
```

For each of the VCF files in the project's S3 Object Store, the dimensions file contains the following metadata:

- VCF filename
- timestamp from when the dimensions of this VCF file was last updated
- the version ID of the latest version of this VCF file in the S3 Object Store
- the number of variants in this VCF file
- the number of samples in this VCF file
- The name of all scaffolds/chromosomes/contigs in this VCF file
- The name of all samples in this VCF file

This information is used by the backend to for the queries. By storing this in a separate VCF dimensions file, the backend does not need to transfer and parse large VCF files every time it needs to fetch any of the above information.

### 2. Query sidecar metadata file for sample IDs and VCF filenames

Attributes of the samples that are typically not part of the data in the VCF file can be stored in a tabular sidecar file where each row describe a sample. The file need to start with a header row that contains one mandatory column: `Sample_ID`. Other columns can be added based on the needs of the project. The mapping between Sample IDs and the VCF file they are contained in is handled by the VCF dimensions file in the backend, so filenames should not be included in the sample metadata file. Since the VCF dimensions file is needed for the queries, the below assumes that the user has run `divbase-cli dimensions update` (as described above) since last time a VCF file was uploaded to the project.

To accomodate queries of user-defined columns, the following query syntax is used:

```bash
divbase-cli query tsv [...] --filter “key1: value1, value2; key2: value3, value4 […]”
```

where `key` is the column name and “value” is the column value to filter on. Values can be stacked by commas, and (intersect) queries across multiple keys can be done by adding a semicolon.

Example: If the test files in `.tests/fixtures/` are uploaded to a project, the command

```bash
divbase-cli query tsv "Area:West of Ireland,Northern Portugal;Sex:F" --metadata-tsv-name sample_metadata.tsv
```

will return the following to the terminal:

```
The results for the query (Area:Northern Portugal):
Unique Sample IDs: ['5a_HOM-I13', '5a_HOM-I14', '5a_HOM-I20', '5a_HOM-I21', '5a_HOM-I7']
Unique filenames: ['HOM_20ind_17SNPs_last_10_samples.vcf.gz', 'HOM_20ind_17SNPs_first_10_samples.vcf.gz']
```

The flag `--show-sample-results` can be added to also display the Sample ID - filename mapping:

```bash
divbase-cli query tsv "Area:Northern Portugal" --metadata-tsv-name sample_metadata.tsv --show-sample-results
```

which, for this example, will return:

```
Name and file for each sample in query results:
Sample ID: '5a_HOM-I7', Filename: 'HOM_20ind_17SNPs_first_10_samples.vcf.gz'
Sample ID: '5a_HOM-I13', Filename: 'HOM_20ind_17SNPs_last_10_samples.vcf.gz'
Sample ID: '5a_HOM-I14', Filename: 'HOM_20ind_17SNPs_last_10_samples.vcf.gz'
Sample ID: '5a_HOM-I20', Filename: 'HOM_20ind_17SNPs_last_10_samples.vcf.gz'
Sample ID: '5a_HOM-I21', Filename: 'HOM_20ind_17SNPs_last_10_samples.vcf.gz'
The results for the query (Area:Northern Portugal):
Unique Sample IDs: ['5a_HOM-I13', '5a_HOM-I14', '5a_HOM-I20', '5a_HOM-I21', '5a_HOM-I7']
Unique filenames: ['HOM_20ind_17SNPs_last_10_samples.vcf.gz', 'HOM_20ind_17SNPs_first_10_samples.vcf.gz']
```

### 3. Query the data in the VCF files based on sample metadata

DivBase uses bcftools to query data contained in the VCF files. One of the main ideas of DivBase is that data in a project can split over multiple VCF files. Queries that require bcftools operations can potentially take long time to complete. Therefore, these queries are submitted to a queue in the DivBase job manager system which is based on Celery. Running `divbase-cli query bcftools-pipe` commands will therefor return a Celery task-id to `stdout`. To check the status of a task:

```bash
divbase-cli query task-status --task-id <ID>
```

If the VCF file query is combined with a sidecar metadata query using `--tsv-filter` which accepts the same the query syntax the sidecar queries described above, DivBase will perform the query on all VCF files identified by the sidecar query, apply any chain bcftools subset/filter operation given by the user, and return a single merged VCF file with the results.

Bcftools subset/filters are set by the `--command` option. The syntax is experimental: it uses existing bcftools commands as defined in the [bcftools manual](https://samtools.github.io/bcftools/bcftools.html) but automatically handles the filenames based on the sidecar metadata query results. At the moment, `view` is the only bcftools command that can be used with `--command`.

At the moment, the sample metadata query needs to be included in the query. To subset on the samples identified in the sample metadata query, use `--command view -s SAMPLES`:

```bash
divbase-cli query bcftools-pipe --tsv-filter "Area:Northern Portugal" --command "view -s SAMPLES" --metadata-tsv-name sample_metadata.tsv
```

When several VCF files are part of the query, the DivBase backend ensure that they each are subset on their own and that all subsets are combined into a single results file using `bcftools merge` or `bcftools concat`, depending on the combinations of samples in each VCF file (see the bcftools manual for more information on the requirements for the `merge` and `concat` commands). Please note that the user do not need to specify `merge` and `concat` or concat in `--command`!

It is possible to create pipes of bcftools operations by semicolon separation. This example first subsets each VCF files on the samples from the metadata queries, and then pipes each results of that subset to a new subset that filters based on scaffold `21` at coordinates `15000000-25000000`:

```bash
divbase-cli query bcftools-pipe --tsv-filter "Area:Northern Portugal" --command "view -s SAMPLES; view -r 21:15000000-25000000" --metadata-tsv-name sample_metadata.tsv
```
